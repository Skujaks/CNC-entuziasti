import math

# Cilindra radius
Rc = 12.5
Hc = 40
# Frezes radius
Rf = 2.5
#Sagataves augstums un radius
Hs = 67
Rs = 25

#mazakie soli uz asim
Cstep = math.radians(1)
Zstep = 0.1

#pozicijas
posX = 0
posZ = 0
posC = 0

#Crot 0 = cw; Crot 1 = ccw
Crot = 1

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
        FrezeSet = ['T060606\n', 'M110\n', 'M146\n', 'G94\n', 'SB=3000 M13\n', 'M'+ str(Crot+15) + '\n', '\n']
        end = ['\n', 'G1 X' + str(round(2 * posX + 20, 3)) + ' Z' +  str(round(posZ, 3)) + '\n', 'G0 35X Z10\n', 'G96S0M3\n', 'M30\n']

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

        #Crot 0 = cw; Crot 1 = ccw
        def CW(C):
            global Crot
            Crot = 0
            file.write('M15\n')

        def CCW():
            global Crot
            Crot = 1
            file.write('M16\n')

        #nomaina 1 uz 0  un  0 uz 1
        def Cflip():
            global Crot
            if Crot == 0:
                Crot = 1

            if Crot == 1:
                Crot = 0
            file.write('M'+ str(Crot+15) + '\n')
        
        def GetA(z):
            a = math.sqrt(Rc**2 + z**2)
            return a
        
        def Theta(a, b):
            #theta = 2*math.atan(b/a)
            R = 2*((b/2)**2 + a**2)
            theta = math.acos((R - b**2)/R)
            return theta
        
        def Gamma(theta):
            gamma = (math.radians(180) - theta)/2
            return gamma
        
        def PosCtoX(a, c):
            r = a*(1/(math.tan(c)))
            if r < 0:
                r = -r
            r += a
            #file.write('B ' + str(round(r, 3)) + '\n')
            return r
        
        #c ir izlabotais c
        def DeltaX(a, c, r):
            if c < 0:
                c = -c
            deltX = r*math.tan(c - math.radians(90))
            if deltX < 0:
                deltX = -deltX
            x2 = PosCtoX(a, c) + deltX
            #file.write('del X ' + str(round(deltX, 3)) + '\n')
            return x2
        
        def DeltaC(a, c, r):
            global posC
            c -= math.radians(90)
            if c < 0:
                c = -c
            x = PosCtoX(a, posC)
            dltC = math.asin(r/x)

            if posC > math.radians(90):
                #file.write('C-\n')
                return -dltC
            else:
                #file.write('C+\n')
                return dltC

        #beta ir grados
        def DeltaZ(beta):
            alpha = math.radians(90 - beta)
            deltZ = Rf*math.cos(alpha)
            return deltZ
        
        # beta ir grados
        def DeltaR(beta):
            alpha = math.radians(90 - beta)
            deltR = Rf*math.sin(alpha)
            return deltR


        def Line(a, b, z, alpha, beta):
            global posX
            global posZ
            global posC

            posZ = z + DeltaZ(beta)
            deltRf = DeltaR(beta)

            theta = Theta(a, b)
            gamma = Gamma(theta)

            for c in range(int(theta/Cstep)+1):
                posC = c*Cstep + gamma
                posC += DeltaC(a, posC, deltRf)
                posX = DeltaX(a, posC, deltRf)
                #file.write('X ' + str(round(posX, 3)) + '\n')
                posC += alpha
                if (posC > 360):
                    posC -= 360

                G101(posX, posZ, posC)

            file.write('\n')
            G1(Rs + 10, posZ)
            Cflip()
            C(alpha + gamma)
            Cflip()
            file.write('\n')

        def Mala1(a, b, Zh, Zoff, alpha, beta):
            for z in range(int(Zh/Zstep)):
                Line(a + z*Zstep*math.tan(math.radians*(beta)), b, -z*Zstep - Zoff, alpha, beta)
#######################################################################################

        file.writelines(setup)

        #Line(Rc, Hc, posZ, 0, 0)
        Mala1()

        file.writelines(end)

if __name__ == '__main__':
    main()