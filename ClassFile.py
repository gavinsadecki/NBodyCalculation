import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#Some of the stuff I think is weird in here. Some of the methods are hacky for the matrix class.
#It assumes certain behavior already that the user has to know and it works here but I would like to fix it and make it more user friendly and general.
#I wanted to get more experience with using objects and classes, which was my intial reason for this rewrite.
#I also felt weird making a polynoial and a matrix class since a polynomial can be represented as an nx1 matrix.
#I started with the polynomial idea since this solution uses a power series approach which is normally represented by polynomails.
#The object G in the calculation was a matrix so I created a matrix object just for that.
class polynomial:
    def __init__(self):
        self.coef = []

    def addTerm(self,nextTerm):
        self.coef.append(nextTerm)
        return self

    def evaluate(self,x):
        result = 0
        for i in range(len(self.coef)):
            result += self.coef[i] * (x**i)
        return result

class planet:
    def __init__(self,name,mass,position,velocity):
        self.name = name
        self.mass = mass
        
        self.posX = position[0]
        self.posY = position[1]
        self.posZ = position[2]

        self.posHistX = [position[0]]
        self.posHistY = [position[1]]
        self.posHistZ = [position[2]]
        
        self.velX = velocity[0]
        self.velY = velocity[1]
        self.velZ = velocity[2]

        self.velHistX = [velocity[0]]
        self.velHistY = [velocity[1]]
        self.velHistZ = [velocity[2]]

    def currentPos(self):
        return [self.posX,self.posY,self.posZ]

    def currentVel(self):
        return [self.velX,self.velY,self.velZ]

    def setPos(self,newPosX,newPosY,newPosZ):
        self.posX = newPosX
        self.posHistX.append(newPosX)
        
        self.posY = newPosY
        self.posHistY.append(newPosY)
                
        self.posZ = newPosZ
        self.posHistZ.append(newPosZ)

    def setVel(self,newVelX,newVelY,newVelZ):
        self.velX = newVelX
        self.velHistX.append(newVelX)
        
        self.velY = newVelY
        self.velHistY.append(newVelY)
                
        self.velZ = newVelZ
        self.velHistZ.append(newVelZ)
        

class matrix:
    def __init__ (self,rows,columns):
        self.values = []
        
        if rows < 1 or columns < 1:
            raise ValueError("Please enter positive non-zero integers only.")
        if type(rows) != int or type(columns) != int:
            raise TypeError("Please enter integers only.")
            
        for i in range(rows):
            tempCol = []
            
            for i in range(columns):
                tempCol.append(0)
                
            self.values.append(tempCol)

    def read(self,row,column):
        try:
            return self.values[row][column]
        except IndexError:
            print('Issue reading value on index provided.')

    def insert(self,row,column,value):
        self.values[row][column].append(value)
        
    def update(self,row,column,value):
        try:
            self.values[row][column] = value
        except IndexError:
            print('Issue updating value on index provided.')

class planetList:
    def __init__(self):
        self.planets = []
        self.gravConst = 6.674*10**-11

    def add(self,planetObj):
        self.planets.append(planetObj)

    def nbody(self,PSI,t,a):
        
        self.G = matrix(len(self.planets),len(self.planets))

        #The power series terms converge to zero really fast. Adding these scale conversions to be able to increase the number, predict the next term and then scale it back down in the end.              
        for p in self.planets:
            p.A = polynomial().addTerm(149597870700*p.posX)
            p.B = polynomial().addTerm(149597870700*p.posY)
            p.C = polynomial().addTerm(149597870700*p.posZ)
            p.D = polynomial().addTerm(1731456.83681*p.velX)
            p.E = polynomial().addTerm(1731456.83681*p.velY)
            p.F = polynomial().addTerm(1731456.83681*p.velZ)

        for indx,p in enumerate(self.planets):
            for indx2,p2 in enumerate(self.planets):
                if indx == indx2:
                    self.G.update(indx,indx2,[0])
                else :
                    Asquared = (p.A.coef[0]-p2.A.coef[0])**2
                    Bsquared = (p.B.coef[0]-p2.B.coef[0])**2
                    Csquared = (p.C.coef[0]-p2.C.coef[0])**2
                    
                    self.G.update(indx,indx2,[(Asquared + Bsquared + Csquared)**(-1/2)])

        for j in range(PSI):
            for indx,p in enumerate(self.planets):
                p.A.addTerm(p.D.coef[j]/(j+1))
                p.B.addTerm(p.E.coef[j]/(j+1))
                p.C.addTerm(p.F.coef[j]/(j+1))

                DTemp = 0
                ETemp = 0
                FTemp = 0
                for indx2,p2 in enumerate(self.planets):
                    GsumTemp = 0
                    
                    for k in range(j+1):
                        for l in range(k+1):
                            for m in range(l+1):
                                for o in range(m+1):
                                    if o <= m and k <= j and l <= k and m <= l:
                                        GsumTemp += self.G.read(indx,indx2)[l-m] * self.G.read(indx,indx2)[m-o] * self.G.read(indx,indx2)[o] * ((p.A.coef[j-k]-p2.A.coef[j-k]) * (p.D.coef[k] - p2.D.coef[k]) + (p.B.coef[j-k] - p2.B.coef[j-k]) * (p.E.coef[k] - p.E.coef[k]) + (p.C.coef[j-k] - p2.C.coef[j-k]) * (p.F.coef[k] - p.F.coef[k]))

                                        base = self.gravConst * p2.mass * self.G.read(indx,indx2)[k-l] * self.G.read(indx,indx2)[l-m] * self.G.read(indx,indx2)[m]
                                                     
                                        DTemp += base * (p2.A.coef[j-k] - p.A.coef[j-k])
                                        ETemp += base * (p2.B.coef[j-k] - p.B.coef[j-k])
                                        FTemp += base * (p2.C.coef[j-k] - p.C.coef[j-k])
                                        
                    self.G.insert(indx,indx2,(-1/(j+1))*GsumTemp)
                                                  
                p.D.addTerm(DTemp/(j+1))
                p.E.addTerm(ETemp/(j+1))
                p.F.addTerm(FTemp/(j+1))
                                                     
            for p in self.planets:
                evalA = p.A.evaluate(t-a)/149597870700
                evalB = p.B.evaluate(t-a)/149597870700
                evalC = p.C.evaluate(t-a)/149597870700
                evalD = p.D.evaluate(t-a)/1731456.83681
                evalE = p.E.evaluate(t-a)/1731456.83681
                evalF = p.F.evaluate(t-a)/1731456.83681

                p.setPos(evalA,evalB,evalC)
                p.setVel(evalD,evalE,evalF)

#End is time past the inital data it will stop, in seconds
#Step is increments in seconds that the estimate will occur for(step = 60 means where the planetary body will be in 60 seconds).
#N is the number of power series terms we will calculate
    def graph(self,end,step,N):
        a = 0
        
        for i in range(0,int(end/step)):
            t = a + step
            self.nbody(N,t,a)
            a += step
        figure = plt.figure()
        actual = figure.add_subplot(111,projection='3d')

        for p in self.planets:
            actual.plot(p.posHistX,p.posHistY,p.posHistZ)
            print(f'Planet Name : {p.name} Xpos : {p.posX}  Ypos : {p.posY}  Zpos : {p.posZ}')
        
        actual.set_xlim(-4,4)
        actual.set_ylim(-4,4)
        actual.set_zlim(-4,4)
        actual.set_xlabel('X')
        actual.set_ylabel('Y')
        actual.set_zlabel('Z')
        plt.show()
