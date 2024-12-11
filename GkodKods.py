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
Cstep = math.radians(2)
Zstep = 0.5

#pozicijas
posX = 0
posZ = 0
posC = 0

maxX = 100

#feed 900
Feed = 900

def main():
    # Specify the name of the output file
    output_file = 'PiemineklisGcode.min'
    # Open the file in write mode ('w')
    with open(output_file, 'w') as file:
        # Write content to the file
        VirpaStart = ['(VIRPA START)\n','\n'
        'G50S3500\n',
        'G96S200M3\n',
        'F0.15\n', '\n']

        VirpaStop = ['(VIRPA STOP)\n','\n'
        'M5\n',
        'G0 X70\n',
        'G00 Z300 X300\n','\n']

        Virpa1 = ['(VIRPA 1)\n','\n'
        'T121212\n',
        'G00 X60 Z0\n'
        'G01 X50\n',
        'X48\n',
        'Z-67\n'
        'X55\n'
        'G00 Z0\n'
        'G01 X47.2\n',
        'Z-47\n',
        'X60\n'
        'G00 Z300 X300\n',
        'T010101\n',
        'X100 Z-47\n'
        'G01 X30\n',
        'X60\n',
        'Z-44\n',
        'X43\n',
        'Z-46\n'
        'X60\n',
        'Z-41\n',
        'X43\n',
        'Z-43\n'
        'X60\n', 
        'Z-38\n',
        'X43\n',
        'Z-40\n'
        'X60\n', 
        'Z-35\n',
        'X43\n',
        'Z-37\n'
        'X60\n', 
        'Z-32\n',
        'X12\n',
        'X60\n', 
        'Z-29\n',
        'X12\n',
        'Z-31\n'
        'X60\n', 
        'Z-28\n',
        'X12\n'
        'Z-30\n'
        'X60\n',
        'G00 Z300 X300\n']

        Virpa2 = ['(VIRPA 2)\n','\n'
        'T010101\n',
        'M76\n',
        'G0 X70 Z-67\n',
        'G96S0M3\n',
        'G1 X0 Z-67\n'
        'X70\n'
        'G0 X300 Z300','\n']

        FrezeStart = ['(FREZE START)\n','\n'
        'T060606\n',
        'M110\n',
        'M146\n',
        'G94\n',
        'SB=3000 M13\n',
        'M16\n',
        'F900\n',
        'G0 X70 Z0\n', '\n']

        FrezeStop = ['(FREZE STOP)\n','\n'
        'M110\n',
        'M146\n',
        'G97\n',
        'SB=0 M13\n',
        'G1 X70\n',
        'G00 Z300 X300\n','\n']

        def maxx(x):
            global maxX
            if x > maxX:
                x = maxX
            return x

        def note(a):
            file.write('(' + str(a) + ')\n' + '\n')

        def G101(x, z, c):
            file.write('G101 X' + str(round(2*maxx(x), 3)) + ' Z' +  str(round(z, 3)) + ' C' + str(round(math.degrees(c), 3)) +'\n')

        def G1(x, z):
            file.write('G1 X' + str(round(2*maxx(x), 3)) + ' Z' +  str(round(z, 3)) +'\n')

        def G0(x, z):
            file.write('G0 X' + str(round(2*maxx(x), 3)) + ' Z' +  str(round(z, 3)) +'\n')
        
        def C(c):
            file.write('C' + str(round(math.degrees(c), 3)) + '\n')
        
        def F(f):
            file.write('F' + str(round(f, 3)) + '\n')

        #Crot 0 = cw; Crot 1 = ccw
        def CW():
            file.write('M15\n')

        def CCW():
            file.write('M16\n')
        
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
            r = a*(1/(math.sin(c) + 0.001))
            if r < 0:
                r = -r
            #r += a
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
            #file.write('A ' + str(round(a, 3)))
            #file.write(' C ' + str(round(c, 3)))
            #file.write(' R ' + str(round(r, 3)) + '\n')
            #dltC = math.asin(r/(x + 0.001))
            dltC = 0
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
            note('Line a' + str(a) + ' b' + str(b) + ' z' + str(z) + ' alpha' + str(alpha) + ' beta' + str(beta))

            global posX
            global posZ
            global posC

            #F(a)

            posZ = z + DeltaZ(beta)
            deltRf = DeltaR(beta)

            theta = Theta(a, b)
            gamma = Gamma(theta)

            for c in range(int(theta/Cstep)+1):
                posC = c*Cstep + gamma
                #C(posC)
                posC += DeltaC(a, posC, deltRf)
                #C(posC)
                posX = DeltaX(a, posC, deltRf)
                #C(posC)
                #file.write('X ' + str(round(posX, 3)) + '\n')
                posC += math.radians(alpha)
                #C(posC)
                if (posC > math.radians(359.9)):
                    posC -= math.radians(360)
                #C(posC)
                #file.write('\n')
                G101(posX, posZ, posC)

            file.write('\n')
            G0(Rs + 10, posZ)
            #CW()
            C(math.radians(alpha + gamma))
            #CCW()
            file.write('\n')

        def Line2(a, b, z, alpha, beta):
            note('Line2 a' + str(a) + ' b' + str(b) + ' Z' + str(z) + ' alpha' + str(alpha) + ' beta' + str(beta))

            global posX
            global posZ
            global posC

            #F(a)

            posZ = z - DeltaZ(beta)
            deltRf = DeltaR(beta)

            theta = Theta(a, b)
            gamma = Gamma(theta)

            Bb = 0
            Bthet = 0
            Bgamma = 0
            Cskip = False

            if a < 6 :
                    Bb = 2*math.sqrt(6**2 - a**2)
                    Bthet = Theta(a, Bb)
                    Bgamma = Gamma(Bthet)

            for c in range(int(theta/Cstep)+1):
                
                posC = c*Cstep + gamma
                #C(posC)
                posC += DeltaC(a, posC, deltRf)


                if posC > gamma - Cstep and posC < Bgamma or posC > Bgamma + Bthet - Cstep and posC < theta + gamma:
                #if posC > gamma and posC < Bgamma or posC > Bgamma + Bthet - Cstep and posC < theta + gamma:
                    #if (posC < math.radians(90)):
                        #note('pirms')
                    #else:
                        #note('pec')
                    #C(posC)
                    posX = DeltaX(a, posC, deltRf)
                    #C(posC)
                    #file.write('X ' + str(round(posX, 3)) + '\n')
                    posC += math.radians(alpha)
                    #C(posC)
                    if (posC > math.radians(359.9)):
                        posC -= math.radians(360)
                    #C(posC)
                    #file.write('\n')
                    G101(posX, posZ, posC)
                    #elif posC > Bgamma + Bthet - Cstep and posC < theta + gamma:
                    #note('pec')
                else:
                    if Cskip == False:
                        #note('else')
                        Cskip = True
                        C(posC + math.radians(alpha) + Cstep)
                        G1(6, posZ)
                        C(Bthet + Bgamma + math.radians(alpha) - 2*Cstep)

            file.write('\n')
            G0(Rs + 10, posZ)
            #CW()
            C(math.radians(alpha + gamma))
            #CCW()
            file.write('\n')


        def Mala1(a, b, Zh, Zoff, alpha, beta):
            note('Mala a' + str(a) + ' b' + str(b) + ' Zh' + str(Zh) + ' Zoff' + str(Zoff) + ' alpha' + str(alpha) + ' beta' + str(beta))
            file.writelines('G0 X70\n')
            G0(35, -Zoff)
            C(math.radians(alpha))
            file.writelines('G1\n')

            for z in range(int(Zh/Zstep)):
                Line(a + z*Zstep*math.tan(math.radians(beta)), b, -z*Zstep - Zoff, alpha, beta)

        def CylTop(R, H, alpha):
            note('Cilindrs top  R' + str(R) + ' H' + str(H) + ' alpha' + str(alpha))
            file.writelines('G0 X70\n')
            G0(35, 0)
            C(math.radians(alpha))
            file.writelines('G1\n')

            for z in range(int(R/Zstep)):
                y = R - z*Zstep
                A = math.acos(y/R)
                Line(R*math.sin(A), H, -(R-y), alpha, math.degrees(math.radians(90) - A))

        def CylB(R, H, alpha):
            note('Cilindrs bottom  R' + str(R) + ' H' + str(H) + ' alpha' + str(alpha))
            file.writelines('G0 X70\n')
            G0(35, -Rc)
            C(math.radians(alpha))
            file.writelines('G1\n')

            for z in range(int(R/Zstep)+1):
                y = z*Zstep
                A = math.acos(y/R)
                Line2(R*math.sin(A), H, -R - y, alpha, math.degrees(math.radians(90) - A))

        def Kaste():
            note('KASTE')
            Mala1(15, 30, 12.5, 32, 0 , 0)
            Mala1(10, 20, 12.5, 32 , 90 , 0)
            Mala1(15, 30, 12.5, 32 , 180 , 0)
            Mala1(10, 20, 12.5, 32 , 270 , 18.435)

        def Cilindrs():
            note('CILINDRS')
            CylTop(Rc, Hc, 90)
            CylTop(Rc, Hc, 180)

            CylB(Rc, Hc, 90)
            CylB(Rc, Hc, 180)

            Mala1(Hc/2, 2*Rc, 2*Rc, 0, 0 , 0)
            Mala1(Hc/2, 2*Rc, 2*Rc, 0, 180, 0)

####################################################################################################
        
        file.writelines('G0 X300 Z300\n\n')

        file.writelines(VirpaStart)
        file.writelines(Virpa1)
        file.writelines(VirpaStop)

        file.writelines(FrezeStart)
        Kaste()
        Cilindrs()
        file.writelines(FrezeStop)

        file.writelines(VirpaStart)
        file.writelines(Virpa2)
        file.writelines(VirpaStop)

        file.writelines('M30')

if __name__ == '__main__':
    main()