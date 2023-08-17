class Patient:
    def __init__(self, name):
        self.name = name
        self.chromosomes = list()

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

        self.HRD = 0
        self.LOH = 0
        self.LST = 0
        self.TAI = 0

    def printAttributes(self):
        print('Patient Name')
        print(self.name)

        print('General CIN')
        print(self.general_cin)
        print('CN-CIN')
        print(self.cn_cin)

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