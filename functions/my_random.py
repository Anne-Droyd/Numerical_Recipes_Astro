import matplotlib.pyplot as plt


class LCG:
    def __init__(self):
        self.instance = 0

    def lcg(self,low = 0, high = 1, seed = None):

        if self.instance == 0:
            if seed is None:
                self.last_output = 42
            else:
                self.last_output = seed
        self.instance +=1
        a =6364136223846793005
        c =9754186451795953191
        m = 2**64
        output = (a*self.last_output+c)%m
        self.last_output = output
        return (high-low)*(output*m**-1)+low



_lcg = LCG()
list = [_lcg.lcg(20,25) for n in range(1000000)]

plt.hist(list,bins=100)
plt.show()