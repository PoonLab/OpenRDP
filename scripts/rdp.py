import numpy as np
import re


class RdpMethod:
    """
    Executes RDP method
    """
    def __init__(self, align, s_names, settings=None, win_size=30, reference=None, min_id=0, max_id=100):
        if settings:
            self.set_options_from_config(settings)
            self.validate_options()

        else:
            self.win_size = win_size
            self.reference = reference
            self.min_id = min_id
            self.max_id = max_id

        self.results = []
        self.s_names = s_names

    def set_options_from_config(self, settings):
        """
        Set the parameters of the RDP method from the config file
        :param settings: a dictionary of settings
        """
        self.win_size = int(settings['window_size'])
        self.reference = settings['reference_sequence']
        self.min_id = int(settings['min_identity'])
        self.max_id = int(settings['max_identity'])

    def validate_options(self):
        """
        Check if the options from the config file are valid
        If the options are invalid, the default value will be used instead
        """
        if self.win_size < 0:
            print("Invalid option for 'window_size'.\nUsing default value (30) instead.")
            self.win_size = 30

        if self.min_id < 0 or self.min_id > 100:
            print("Invalid option for 'min_identity'.\nUsing default value (0) instead.")
            self.min_id = 0

        if self.max_id < 0 or self.max_id > 100:
            print("Invalid option for 'max_identity'.\nUsing default value (100) instead.")
            self.min_id = 100

    @staticmethod
    def remove_uninformative_sites(triplet):
        """
        Remove sites that are all the same or all different
        :param triplet: n x 3 numpy array where n is the length of the alignment
        :return: a tuple containing the informative sites and the positions of informative sites
        """
        infor_sites = []
        uninfor_sites = []
        # Find positions of sites that are all the same sites or all sites that are different
        for i in range(triplet.shape[1]):
            col = triplet[:, i]
            if np.unique(col).shape[0] == 2:
                infor_sites.append(i)
            else:
                uninfor_sites.append(i)

        # Build "new alignment"
        new_trp = triplet[:, infor_sites]

        return new_trp, infor_sites, uninfor_sites

    def execute(self, alignment, triplet, G):
        """
        Performs RDP detection method for one triplet of sequences
        :param alignment: the Alignment object
        :param triplet: the triplet
        :param G: the number of possible triplets
        :return: the coordinates of the potential recombinant region and the p_value
        """
        # Get the triplet sequences
        trp_seqs = []
        for seq_num in triplet:
            trp_seqs.append(alignment[seq_num])
        trp_seqs = np.array(trp_seqs)

        rec_name = self.s_names[triplet[0]]
        p1_name = self.s_names[triplet[1]]
        p2_name = self.s_names[triplet[2]]

        # Find the informative sites
        new_trp, sites, unsites = self.remove_uninformative_sites(trp_seqs)

        # Get the three pairs of sequences
        ab = np.array([new_trp[0], new_trp[1]])
        bc = np.array([new_trp[1], new_trp[2]])
        ac = np.array([new_trp[0], new_trp[2]])

        len_trp = new_trp.shape[1]

        # 2. Sliding window over subsequence and calculate average percent identity at each position
        recombinant_regions = ''  # Recombinant regions denoted by ones
        coord = []
        for i in range(len_trp - self.win_size):
            reg_ab = ab[:, i: self.win_size + i]
            reg_bc = bc[:, i: self.win_size + i]
            reg_ac = ac[:, i: self.win_size + i]

            # Calculate percent identity in each window
            a_b, b_c, a_c = 0, 0, 0
            for j in range(reg_ab.shape[1]):
                if reg_ab[0, j] == reg_ab[1, j]:
                    a_b += 1
                if reg_bc[0, j] == reg_bc[1, j]:
                    b_c += 1
                if reg_ac[0, j] == reg_ac[1, j]:
                    a_c += 1

            percent_identity_ab = a_b / len_trp * 100
            percent_identity_bc = b_c / len_trp * 100
            percent_identity_ac = a_c / len_trp * 100

            # Identify recombinant regions
            if percent_identity_ac > percent_identity_ab or percent_identity_bc > percent_identity_ab:
                recombinant_regions += "1"
                coord.append(i)
            else:
                recombinant_regions += "0"

        # 3. Record significance of events
        recomb_idx = [(m.span()) for m in re.finditer('1+', recombinant_regions)]

        # Convert coordinates from  window-level to alignment-level and record number of windows
        coords = []
        for x, y in recomb_idx:
            coords.append((sites[x], sites[y - 1]))

        for coord in coords:
            n = coord[1] - coord[0]     # Length of putative recombinant region

            if n > 0:
                # m is the proportion of nts in common between either A or B and C in the recombinant region
                nts_in_a = trp_seqs[0][coord[0]: coord[1]]
                nts_in_c = trp_seqs[2][coord[0]: coord[1]]
                m = 0
                for i in range(n):
                    if nts_in_a[i] == nts_in_c[i]:
                        m += 1

                # p is the proportion of nts in common between either A or B and C in the entire subsequence
                id_in_seq = 0
                for j in range(trp_seqs.shape[1]):
                    if trp_seqs[0][j] == trp_seqs[2][j]:
                        id_in_seq += 1
                p = id_in_seq / trp_seqs.shape[1]

                # Calculate p_value
                val = 0
                log_n_fact = np.sum(np.log(np.arange(1, n+1)))  # Convert to log space to prevent integer overflow
                for i in range(m, n):
                    log_i_fact = np.sum(np.log(np.arange(1, i+1)))
                    log_ni_fact = np.sum(np.log(np.arange(1, n-i+1)))
                    try:
                        val += np.math.exp((log_n_fact - (log_i_fact + log_ni_fact)) + np.log(p**n) + np.log((1-p)**(n-i)))
                    except ZeroDivisionError:
                        pass

                uncorr_pvalue = (len_trp / n) * val
                corr_p_value = G * uncorr_pvalue

            else:
                uncorr_pvalue = 'NS'
                corr_p_value = 'NS'

            self.results.append((rec_name, p1_name, p2_name, coord, uncorr_pvalue, corr_p_value))

        return
