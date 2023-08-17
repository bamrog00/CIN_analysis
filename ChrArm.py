class ChrArm:
    def __init__(self, arm, chromosome_number, start, end, length):
        self.chromosome_number = chromosome_number
        self.chromosome_arm = chromosome_number + '_' + arm
        self.start = start
        self.end = end
        self.length = length

        self.status = 'None'
        self.average_cn = 0
        self.general_cin = 0
        self.cn_cin = 0

        self.number_cn_amp = 0
        self.number_cn_gain = 0
        self.number_cn_loh = 0
        self.number_homo_del = 0
        self.number_hemi_del = 0
        self.number_neutral = 0
        self.number_gain = 0
        self.number_amp = 0
        self.mean_seg_length = 0
        self.mean_cn_loh_length = 0
        self.mean_cn_gain_length = 0
        self.mean_cn_amp_length = 0
        self.mean_homo_del_length = 0
        self.mean_hemi_del_length = 0
        self.mean_neutral_length = 0
        self.mean_gain_length = 0
        self.mean_amp_length = 0

        self.segments_indices = list()
        self.segments_coverage = list()
        self.segments_length = list()
        self.segments_event = list()
        self.segments_copy_number = list()

    def printAttributes(self):
        print('Chr number')
        print(self.chromosome_number)
        print('Chromosome arm combined name')
        print(self.chromosome_arm)
        print('Start')
        print(self.start)
        print('End')
        print(self.end)
        print('Lenght')
        print(self.length)

        print('Status')
        print(self.status)
        print('General CIN')
        print(self.general_cin)
        print('CN-CIN')
        print(self.cn_cin)
        print('Average Copy number')
        print(self.average_cn)

        print('Number CN-LOH')
        print(self.number_cn_loh)
        print('Number CN-GAIN')
        print(self.number_cn_gain)
        print('Number CN-AMP')
        print(self.number_cn_amp)
        print('Number Homo-Del')
        print(self.number_homo_del)
        print('Number Hemi-Del')
        print(self.number_hemi_del)
        print('Number Neutral')
        print(self.number_neutral)
        print('Number Gain')
        print(self.number_gain)
        print('Number Amp')
        print(self.number_amp)
        print('Mean segment length')
        print(self.mean_seg_length)
        print('Mean CN-LOH length')
        print(self.mean_cn_loh_length)
        print('Mean CN-GAIN length')
        print(self.mean_gain_length)        
        print('Mean CN-AMP length')
        print(self.mean_amp_length)
        print('Mean Homo-Del length')
        print(self.mean_homo_del_length)
        print('Mean Hemi-Del length')
        print(self.mean_hemi_del_length)
        print('Mean Neutral length')
        print(self.mean_neutral_length)
        print('Mean Gain length')
        print(self.mean_gain_length)
        print('Mean Amp length')
        print(self.mean_amp_length)

        print('Indices Segments')
        print(self.segments_indices)
        print('Coverage Segments')
        print(self.segments_coverage)
        print('Length segments')
        print(self.segments_length)
        print('Event segements')
        print(self.segments_event)
        print('Segments copy number')
        print(self.segments_copy_number)