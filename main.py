import matplotlib.pyplot as plt
import sys
import math
import numpy as np


def read_xyz_file(path, first_line_offset=2):
    file = open(path, "r")
    slines = file.readlines()[first_line_offset:]
    return [[line.split()[0] for line in slines],
            [[float(coord) for coord in line.split()[1:]] for line in slines]]


def dist(v1, v2):
    return math.sqrt(math.pow(v1[0] - v2[0], 2) + math.pow(v1[1] - v2[1], 2) + math.pow(v1[2] - v2[2], 2))


#   Возвращает центр масс [x,y,z].
#   Входные параметры: массив координат [][x,y,z] и массив масс []
def calc_center_of_mass(coords, masses):
    cm = [0, 0, 0]
    for i in range(len(coords)):
        cm[0] += coords[i][0] * masses[i]
        cm[1] += coords[i][1] * masses[i]
        cm[2] += coords[i][2] * masses[i]
    return [c / sum(masses) for c in cm]


#   Возвращает массив расстояний |r_cm - r_i|.
#   Входные параметры: массив координат [][x,y,z] и массив масс []
def calc_d_array(coords, masses):
    cm = calc_center_of_mass(coords, masses)
    return sorted([dist(coord, cm) for coord in coords])


def gauss(x, sigma, mu):
    return np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (sigma * np.sqrt(2.0 * np.pi))


#   Строит массив точек x, y для отображения [[],[]]
#   Входные параметры: массив расстояний [], левая и правая граница отображения по x, величина шага для отображения,
#   сигма для гауссовской функции
def calc_xy_arrays(d_array, r_min, r_max, step, sigma):
    r_len = r_max - r_min
    npoints = int(r_len / step) + 1
    x_arr = np.linspace(r_min, r_max, npoints)
    y_arr = [0] * npoints
    for i in range(npoints):
        for pt in d_array:
            y_arr[i] += gauss(x_arr[i], sigma, pt)
    return [x_arr, y_arr]


def find_r_min_max(d_arrays):
    g_min = 0
    g_max = 0
    for d_array in d_arrays:
        g_min = min(g_min, min(d_array))
        g_max = max(g_max, max(d_array))
    return [g_min, g_max]


#   Вычисляет массив расстояний для каждого заданного в аргументах файла [][].
def calc_d_arrays(file_names):
    d_arrays = []
    for filename in file_names:
        read_result = read_xyz_file(filename)
        coords = read_result[1]  # координаты частиц [][x,y,z]
        masses = [masses_lib[mass_id] for mass_id in read_result[0]]  # массы частиц []
        d_array = calc_d_array(coords, masses)  # массив расстояний от центра массы до частиц |r_cm - r_i|
        d_arrays.append(d_array)
    return d_arrays


########################################################################################################################

masses_lib = {
    'Ir': 192.217,
    'Ag': 107.8682,
    'Au': 196.966
}

# file_names = ["example1.xyz", "example2.xyz"]
file_names = ["example2.xyz", "example3.xyz", "example4.xyz"]

if __name__ == "__main__":
    if len(sys.argv) > 1:
        file_names = sys.argv[1:]
    # else:
    # raise SystemExit

d_arrays = calc_d_arrays(file_names)
r_min_max = find_r_min_max(d_arrays)
xy_result_arrays = [calc_xy_arrays(d_array, r_min_max[0], r_min_max[1], 0.01, 0.01) for d_array in d_arrays]

for xy_ra in xy_result_arrays:
    plt.plot(xy_ra[0], xy_ra[1], '-')

plt.legend(file_names)
plt.show()
