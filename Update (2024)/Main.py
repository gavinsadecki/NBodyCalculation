from Planets import Planet, Planets


def f(x):
    return 10**x


# Data comes from JPL horizons.
# These represent the planets at a set date and time that I don't recall since I grabbed the data from 2018.
sun = Planet('Sun',
             1.988544*f(30),
             [0, 0, 0],
             [0, 0, 0])
mercury = Planet('Mercury',
                 3.302*f(23),
                 [-1.407280799640445*f(-1), -4.439009577663330*f(-1), -2.334555971312334*f(-2)],
                 [2.116887137167173*f(-2), -7.097975438870807*f(-3), -2.522830951443754*f(-3)])
venus = Planet('Venus',
               4.8685*f(24),
               [-7.186302169039649*f(-1), -2.250380105571625*f(-2), 4.117184137682463*f(-2)],
               [5.135327579269579*f(-4), -2.030614162239802*f(-2), -3.071745100210852*f(-4)])
earth = Planet('Earth',
               5.97219*f(24),
               [-1.685246489174995*f(-1), 9.687833048228511*f(-1), -4.120973411130758*f(-6)],
               [-1.723394583068879*f(-2), -3.007660259271771*f(-3), 3.562931614781975*f(-8)])
mars = Planet('Mars',
              6.4185*f(23),
              [1.390361066087240, -2.100972123734226*f(-2), -3.461801385164819*f(-2)],
              [7.479271243289054*f(-4), 1.518629867736057*f(-2), 2.997531995727463*f(-4)])
jupiter = Planet('Jupiter',
                 1.89813*f(27),
                 [4.003460455020136, 2.935353231223561, -1.018232818649105*f(-1)],
                 [-4.563750862068905*f(-3), 6.447274190510942*f(-3), 7.546991834715161*f(-5)])
saturn = Planet('Saturn',
                5.68319*f(26),
                [6.408556034969109, 6.568042753804105, -3.691272981086989*f(-1)],
                [-4.290540499096270*f(-3), 3.891990892070301*f(-3), 1.026097570670085*f(-4)])
uranus = Planet('Uranus',
                8.68103*f(25),
                [1.443051758274212*f(1), -1.373565828163257*f(1), -2.381283040241605*f(-1)],
                [2.678466083136532*f(-3), 2.672427506613704*f(-3), -2.474809404498369*f(-5)])
neptune = Planet('Neptune',
                 1.0241*f(26),
                 [1.681075678673751*f(1), -2.499265126481678*f(1), 1.272710177158801*f(-1)],
                 [2.579215311143540*f(-3), 1.776354511745528*f(-3), -9.619957499965839*f(-5)])

planets = Planets([sun, mercury, venus, earth, mars, jupiter, saturn, uranus, neptune])
planets.graph(31557600, 3600, 3)