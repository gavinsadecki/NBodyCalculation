import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


class Planet:
    def __init__(self, name: str, mass: float, position: list, velocity: list):
        """
        :param name (str): Name of the planetary body, used for labeling and tracking.
        :param mass (float): Mass of the body in kilograms.
        :param position (list(float)): x,y,z coordinates of the body. Units are in astronomical units(au).
        :param velocity (list(float)): x,y,z velocities of the body. Units are in au/d(ay).
        """
        if len(position) != 3:
            raise ValueError("Position must be 3 dimensional")
        if len(velocity) != 3:
            raise ValueError("Velocity must be 3 dimensional")

        self.name = name
        self.mass = mass

        self.pos_x, self.pos_y, self.pos_z = position
        self.vel_x, self.vel_y, self.vel_z = velocity

        self.pos_hist_x, self.pos_hist_y, self.pos_hist_z = [[position[0]], [position[1]], [position[2]]]
        self.vel_hist_x, self.vel_hist_y, self.vel_hist_z = [[velocity[0]], [velocity[1]], [velocity[2]]]

    def get_position(self):
        return [self.pos_x, self.pos_y, self.pos_z]

    def get_velocity(self):
        return [self.vel_x, self.vel_y, self.vel_z]

    def set_position(self, new_pos_x: float, new_pos_y: float, new_pos_z: float):
        self.pos_x, self.pos_y, self.pos_z = new_pos_x, new_pos_y, new_pos_z

        self.pos_hist_x.append(new_pos_x)
        self.pos_hist_y.append(new_pos_y)
        self.pos_hist_z.append(new_pos_z)

    def set_velocity(self, new_vel_x: float, new_vel_y: float, new_vel_z: float):
        self.vel_x, self.vel_y, self.vel_z = new_vel_x, new_vel_y, new_vel_z

        self.vel_hist_x.append(new_vel_x)
        self.vel_hist_y.append(new_vel_y)
        self.vel_hist_z.append(new_vel_z)


class Planets:
    def __init__(self, planets: list = None):
        if planets is None:
            self.list = []
        else:
            self.list = planets

    def add(self, planet: Planet):
        self.list.append(planet)

    def nbody(self, terms: int, t: int, a: int):
        const_position = 149597870700
        const_velocity = 1731456.83681
        const_gravity = 6.674 * 10 ** -11

        term_g = np.ndarray((len(self.list), len(self.list), terms + 1))
        term_g.fill(0)

        # The power series terms converge to zero really fast.
        # Adding these scale conversions to be able to increase the number,
        # predict the next term and then scale it back down in the end.
        for p in self.list:
            p.terms_a = [const_position * p.pos_x]
            p.terms_b = [const_position * p.pos_y]
            p.terms_c = [const_position * p.pos_z]
            p.terms_d = [const_velocity * p.vel_x]
            p.terms_e = [const_velocity * p.vel_y]
            p.terms_f = [const_velocity * p.vel_z]

        for index, p in enumerate(self.list):
            for index2, p2 in enumerate(self.list):
                if index != index2:
                    squared_a = (p.terms_a[0] - p2.terms_a[0]) ** 2
                    squared_b = (p.terms_b[0] - p2.terms_b[0]) ** 2
                    squared_c = (p.terms_c[0] - p2.terms_c[0]) ** 2

                    term_g[index][index2][0] = np.pow(squared_a + squared_b + squared_c, -0.5)

        for j in range(terms):
            for index, p in enumerate(self.list):
                p.terms_a.append(p.terms_d[j] / (j + 1))
                p.terms_b.append(p.terms_e[j] / (j + 1))
                p.terms_c.append(p.terms_f[j] / (j + 1))

                temp_d = 0
                temp_e = 0
                temp_f = 0
                for index2, p2 in enumerate(self.list):
                    temp_g = 0

                    for k in range(j + 1):
                        for l in range(k + 1):
                            for m in range(l + 1):
                                for o in range(m + 1):
                                    if o <= m <= l <= k <= j:
                                        temp_g += (
                                                        term_g[index, index2, l-m] *
                                                        term_g[index, index2, m-o] *
                                                        term_g[index, index2, o] *
                                                        (
                                                            (p.terms_a[j - k] - p2.terms_a[j - k]) *
                                                            (p.terms_d[k] - p2.terms_d[k]) +
                                                            (p.terms_b[j - k] - p2.terms_b[j - k]) *
                                                            (p.terms_e[k] - p.terms_e[k]) +
                                                            (p.terms_c[j - k] - p2.terms_c[j - k]) *
                                                            (p.terms_f[k] - p.terms_f[k])
                                                        )
                                                    )

                                        base = const_gravity * p2.mass * term_g[index, index2, k - l] * \
                                            term_g[index][index2][l-m] * term_g[index, index2, m]

                                        temp_d += base * (p2.terms_a[j - k] - p.terms_a[j - k])
                                        temp_e += base * (p2.terms_b[j - k] - p.terms_b[j - k])
                                        temp_f += base * (p2.terms_c[j - k] - p.terms_c[j - k])

                    term_g[index, index2, j + 1] = (-1 / (j + 1)) * temp_g

                p.terms_d.append(temp_d / (j + 1))
                p.terms_e.append(temp_e / (j + 1))
                p.terms_f.append(temp_f / (j + 1))

            for p in self.list:
                eval_a = calculate_reverse_polynomial(p.terms_a, t-a) / const_position
                eval_b = calculate_reverse_polynomial(p.terms_b, t-a) / const_position
                eval_c = calculate_reverse_polynomial(p.terms_c, t-a) / const_position
                eval_d = calculate_reverse_polynomial(p.terms_d, t-a) / const_velocity
                eval_e = calculate_reverse_polynomial(p.terms_e, t-a) / const_velocity
                eval_f = calculate_reverse_polynomial(p.terms_f, t-a) / const_velocity

                p.set_position(eval_a, eval_b, eval_c)
                p.set_velocity(eval_d, eval_e, eval_f)

    def graph(self, end: int, step: int, terms: int):
        """
        :param end: End time to model the planetary movement to
        :param step: Steps taken per calculation. Smaller steps improve accuracy but increase computation
        :param terms: Number of power series terms calculated, more terms means more accuracy
        :return: None
        """
        for a in range(int(end / step)):
            self.nbody(terms, (a + 1) * step, a * step)

        figure = plt.figure()
        graph = figure.add_subplot(111, projection='3d')

        for planet in self.list:
            graph.plot(planet.pos_hist_x, planet.pos_hist_y, planet.pos_hist_z)

        graph.set_xlim(-4, 4)
        graph.set_ylim(-4, 4)
        graph.set_zlim(-4, 4)
        graph.set_xlabel('X')
        graph.set_ylabel('Y')
        graph.set_zlabel('Z')
        plt.show()


def calculate_reverse_polynomial(terms: list, variable: float):
    """
    Numpy polynomial evaluate uses the last index as the first term of the polynomial
    Terms are inserted with the first term as the first value in the list, so using this method
    to allow numpy to calculate the polynomial
    :param terms: List of polynomial terms
    :param variable: Variable to use in the calculation of the polynomial
    :return: Float
    """
    terms.reverse()
    result = np.polyval(terms, variable)
    terms.reverse()

    return result