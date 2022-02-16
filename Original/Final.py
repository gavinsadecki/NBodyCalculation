import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

def nbody(P,V,M,PSI,t,a):
    A,B,C,D,E,F,G,GC,n = [],[],[],[],[],[],[],(6.674*10**-11),len(P)
    sumAL,sumBL,sumCL,sumDL,sumEL,sumFL = [],[],[],[],[],[]
    for i in range(0,n):
        A.append([149597870700*P[i][0]])
        B.append([149597870700*P[i][1]])
        C.append([149597870700*P[i][2]])
        D.append([1731456.83681*V[i][0]])
        E.append([1731456.83681*V[i][1]])
        F.append([1731456.83681*V[i][2]])
    for i in range(0,n):
        for j in range(0,n):
            if i == j :
                G.append([0])
            else:
                G.append([((A[i][0]-A[j][0])**2+(B[i][0]-B[j][0])**2+(C[i][0]-C[j][0])**2)**(-1/2)])
    for j in range(0,PSI):        
        for r in range(0,n):
            A[r].append(D[r][j]/(j+1))
            B[r].append(E[r][j]/(j+1))
            C[r].append(F[r][j]/(j+1))           
        for i in range(0,len(G)):
            GTemp,first = 0,int(i/n)
            second = i-first*n
            for k in range(0,j+1):
                for l in range(0,k+1):
                    for m in range(0,l+1):
                        for o in range(0,m+1):
                            if o<m+1 and k<j+1 and l<k+1 and m<l+1:
                                GTemp += ((G[i][l-m]*G[i][m-o]*G[i][o])*(A[first][j-k]*D[first][k]-A[first][j-k]*D[second][k]-A[second][j-k]*D[first][k]+A[second][j-k]*D[second][k]+B[first][j-k]*E[first][k]-B[first][j-k]*E[second][k]-B[second][j-k]*E[first][k]+B[second][j-k]*E[second][k]+C[first][j-k]*F[first][k]-C[first][j-k]*F[second][k]-C[second][j-k]*F[first][k]+C[second][j-k]*F[second][k]))
            G[i].append((-1/(j+1))*GTemp)      
        for i in range(0,n):
            DTemp,ETemp,FTemp=0,0,0
            for r in range(0,n):
                for k in range(0,j+1):
                    for l in range(0,k+1):
                        for m in range(0,l+1):
                            for o in range(0,m+1):
                                if o<m+1 and k<j+1 and l<k+1 and m<l+1:
                                    if i!=r:
                                        DTemp,ETemp,FTemp = DTemp + ((GC*M[r]*G[i*n+r][k-l]*G[i*n+r][l-m]*G[i*n+r][m])*(A[r][j-k]-A[i][j-k])),ETemp + ((GC)*(M[r]*G[i*n+r][k-l]*G[i*n+r][l-m]*G[i*n+r][m])*(B[r][j-k]-B[i][j-k])),FTemp + ((GC)*(M[r]*G[i*n+r][k-l]*G[i*n+r][l-m]*G[i*n+r][m])*(C[r][j-k]-C[i][j-k]))                  
            D[i].append((1/(j+1))*DTemp)
            E[i].append((1/(j+1))*ETemp)
            F[i].append((1/(j+1))*FTemp)
    for i in range(0,n):
        sumA,sumB,sumC,sumD,sumE,sumF = 0,0,0,0,0,0
        for j in range(0,PSI):
            sumA,sumB,sumC,sumD,sumE,sumF = sumA + A[i][j]*((t-a)**j),sumB + B[i][j]*((t-a)**j),sumC + C[i][j]*((t-a)**j),sumD + D[i][j]*((t-a)**j),sumE + E[i][j]*((t-a)**j), sumF + F[i][j]*((t-a)**j)
        sumAL.append(sumA/149597870700)
        sumBL.append(sumB/149597870700)
        sumCL.append(sumC/149597870700)
        sumDL.append(sumD/1731456.83681)
        sumEL.append(sumE/1731456.83681)
        sumFL.append(sumF/1731456.83681)       
    FinalPos,FinalVel = [],[]
    for i in range(0,n):
        FinalPos.append([sumAL[i],sumBL[i],sumCL[i]])
        FinalVel.append([sumDL[i],sumEL[i],sumFL[i]])
    return FinalPos,FinalVel

def f(x):
    return 10**x


#[Sun,Mercury,Venus,Earth,Mars,Jupiter,Saturn,Uranus,Neptune,Moon]
M = [1.988544*f(30),3.302*f(23),48.685*f(23),5.97219*f(24),6.4185*f(23),1898.13*f(24),5.68319*f(26),86.8103*f(24),102.41*f(24)]
P = [[0,0,0],[-1.407280799640445*f(-1),-4.439009577663330*f(-1),-2.334555971312334*f(-2)],[-7.186302169039649*f(-1),-2.250380105571625*f(-2),4.117184137682463*f(-2)],[-1.685246489174995*f(-1),9.687833048228511*f(-1),-4.120973411130758*f(-6)],[1.390361066087240,-2.100972123734226*f(-2),-3.461801385164819*f(-2)],[4.003460455020136,2.935353231223561,-1.018232818649105*f(-1)],[6.408556034969109,6.568042753804105,-3.691272981086989*f(-1)],[1.443051758274212*f(1),-1.373565828163257*f(1),-2.381283040241605*f(-1)],[1.681075678673751*f(1),-2.499265126481678*f(1),1.272710177158801*f(-1)]] #[-1.706480096344460*f(-1),9.671664156658847*f(-1),2.402396170914541*f(-4)]
V = [[0,0,0],[2.116887137167173*f(-2),-7.097975438870807*f(-3),-2.522830951443754*f(-3)],[5.135327579269579*f(-4),-2.030614162239802*f(-2),-3.071745100210852*f(-4)],[-1.723394583068879*f(-2),-3.007660259271771*f(-3),3.562931614781975*f(-8)],[7.479271243289054*f(-4),1.518629867736057*f(-2),2.997531995727463*f(-4)],[-4.563750862068905*f(-3),6.447274190510942*f(-3),7.546991834715161*f(-5)],[-4.290540499096270*f(-3),3.891990892070301*f(-3),1.026097570670085*f(-4)],[2.678466083136532*f(-3),2.672427506613704*f(-3),-2.474809404498369*f(-5)],[2.579215311143540*f(-3),1.776354511745528*f(-3),-9.619957499965839*f(-5)]] #[-1.691006420027502*f(-2),-3.469659799257778*f(-3),-9.111637218497984*f(-7)]

def Error(Z,Actual):
    return max([(Actual[0]-Z[0][1][0])/Actual[0],(Actual[1]-Z[0][1][1])/Actual[1],(Actual[2]-Z[0][1][2])/Actual[2]])


#end = end time point in seconds, step is the step between each interval
def graphing(P,V,M,end,step,N):
    T1 = time.time()
    a=0
    x,y,z = [],[],[]

    for i in range(0,len(P)):
        x.append([P[i][0]])
        y.append([P[i][1]])
        z.append([P[i][2]])
    
    for i in range(0,int(end/step)):
        t = a+step
        Z = nbody(P,V,M,N,t,a)
        a += step

        for j in range(0,len(P)):
            P[j][0] = Z[0][j][0]
            P[j][1] = Z[0][j][1]
            P[j][2] = Z[0][j][2]
            V[j][0] = Z[1][j][0]
            V[j][1] = Z[1][j][1]
            V[j][2] = Z[1][j][2]
            x[j].append(Z[0][j][0])
            y[j].append(Z[0][j][1])
            z[j].append(Z[0][j][2])
    T2 = time.time()
    print(T2 - T1)
    figure = plt.figure()
    actual = figure.add_subplot(111,projection='3d')
    Color = ['k','r','m','g','r','y','k','b','c']
    for i in range(0,len(P)):
        actual.plot(x[i],y[i],z[i],Color[i])
        print(f' Xpos : {x[i][-1]}  Ypos : {y[i][-1]}  Zpos : {z[i][-1]}')

    
    actual.set_xlim(-4,4)
    actual.set_ylim(-4,4)
    actual.set_zlim(-4,4)
    actual.set_xlabel('X')
    actual.set_ylabel('Y')
    actual.set_zlabel('Z')
    plt.show()

graphing(P,V,M,1*31104000,21600,3)
