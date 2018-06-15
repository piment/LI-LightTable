from modules import opc

class LightTable:

    def __init__(self):
        self.client = opc.Client('localhost:7890')
        self.num_LEDs = 512

        self.AllOff = [(0, 0, 0)] * self.num_LEDs
        self.AllOn = [(255, 255, 255)] * self.num_LEDs
        self.Error = [(255, 0, 0)] * self.num_LEDs
        self.lines = [(0, 0, 0)] * 512

        self.color_red = (255, 0, 0)
        self.color_blue = (0, 0, 255)
        self.color_green = (0, 255, 100)
        self.color_white = (255, 255, 255)
        self.color_off = (0, 0, 0)

        self.pattern_case_4 = [
            [(8, 9, 72, 73), (6, 7, 70, 71), (4, 5, 68, 69), (2, 3, 66, 67), (0, 1, 64, 65)],
            [(136, 137, 200, 201), (134, 135, 198, 199), (132, 133, 196, 197), (130, 131, 194, 195),
             (128, 129, 192, 193)],
            [(264, 265, 328, 329), (262, 263, 326, 327), (260, 261, 324, 325), (258, 259, 322, 323),
             (256, 257, 320, 321)],
            [(392, 393, 456, 457), (390, 391, 454, 455), (388, 389, 452, 453), (386, 387, 450, 451),
             (384, 385, 448, 449)]
        ]

    def get_case(self, num):
        if num == 0 or num > 20:
            return (0,0,0,0)

        if num <= 5:  # num = 1 a 5
            return self.pattern_case_4[0][num - 1]

        if 5 < num <= 10:
            return self.pattern_case_4[1][num - 6]

        if 10 < num <= 15:
            return self.pattern_case_4[2][num - 11]

        if 15 < num <= 20:
            return self.pattern_case_4[3][num - 16]

    def table_raz(self):
        self.lines = self.lines = [(0, 0, 0)] * 512

    def set_case(self, case, color):
        for pos in self.get_case(case):
            self.lines[pos] = color
