import math

# Cilindra radius
Rc = 12.5
Hc = 40
# Frezes radius
Rd = 2.5
#Sagataves augstums un radius
Hs = 67
Rs = 25

Cstep = math.radians(1)
Zstep = 0.1

posX = 0
posZ = 0
posC = 0

#feed 900
Feed = 900


def main():
    # Specify the name of the output file
    output_file = 'PiemineklisGcode.min'
    # Open the file in write mode ('w')
    with open(output_file, 'w') as file:
        # Write content to the file
        setup = ['T121212\n', 'G50S3500\n', 'G96S200M3\n', 'G0\n', '\n']
        VirpaSet = ['T121212\n', 'G50S3500\n', 'G96S200M3\n', '\n']
        VirpaStop = ['T121212\n', 'G50S3500\n', 'G96S200M3\n', '\n']
        Virpa1 = ['G1 X50 G1 Z-67\n' 'G1 X48\n' 'G1 Z-0\n' 'G1 X47.2\n' 'G1 Z-47\n' 'G1 X43\n' 'G1 Z-25\n' 'G1 X12\n' 'G1 Z-32\n', '\n']
        FrezeSet = ['T060606\n', 'M110\n', 'M146\n', 'G94\n', 'SB=3000 M13\n', 'M15\n', '\n']
        end = ['\n', 'G1 X' + str(round(2*posX + 20, 3)) + ' Z' +  str(round(posZ, 3)) + '\n', 'G0 35X Z10\n', 'G96S0M3\n', 'M30\n']

        def G101(x, z, c):
            file.write('G101 X' + str(round(2*x, 3)) + ' Z' +  str(round(z, 3)) + ' C' + str(round(math.degrees(c), 3)) +'\n')

        def G1(x, z):
            file.write('G1 X' + str(round(2*x, 3)) + ' Z' +  str(round(z, 3)) +'\n')

        def G0(x, z):
            file.write('G0 X' + str(round(2*x, 3)) + ' Z' +  str(round(z, 3)) +'\n')
        
        def C(c):
            file.write('C' + str(round(math.degrees(c), 3)) + '\n')
        
        def F(f):
            file.write('F' + str(round(f, 3)) + '\n')

        def CW():
            file.write('M15\n')

        def CCW():
            file.write('M16\n')

        def Beta(z):
            beta = math.asin(z/Rc/2)
            return beta
        
        def GetA(z):
            a = math.sqrt(Rc**2 + z**2)
            return a
        
        def Theta(a, b):
            #theta = 2*math.atan(b/a)
            Rmax = 2*((b/2)**2 + a**2)
            theta = math.acos((Rmax - b**2)/Rmax)
            return theta
        
        def Gamma(theta):
            gamma = (math.radians(180) - theta)/2
            return gamma
        
        def Pos():
            return 0
        
        def DeltaX(a, b, z, c, case):
            x = a
            return x
        
        def DeltaZ(a, b):
            z = b
            return z

        def LineA(a, b, z, alpha, case):
            posZ = z
            theta = Theta(a, b)
            gamma = Gamma(theta)

            for c in range(int(theta/Cstep)+2):
                r = a*(1/(math.tan(c*Cstep+gamma)))
                if r < 0:
                    r = -r
                r += a
                if (c*Cstep + alpha + gamma > 360):
                    posC = c*Cstep + alpha + gamma - 360
                else:
                    posC = c*Cstep + alpha + gamma
                #case 0 = taisna virsma; 1 = slipa virsma; 2 = cilindrs
                posX = DeltaX(r, z, posC, posZ, case)
                posZ = DeltaZ(r, z)
                G101(posX, posZ, posC)
            G1(Rs + 10, posZ)
            C(alpha + gamma)
           

#######################################################################################

        G1(30, 4)
        C(70)

        file.writelines(setup)

        LineA(Rc, Hc, posZ, 0, 0)
        #CylTop()
        file.writelines(end)

if __name__ == '__main__':
    main()