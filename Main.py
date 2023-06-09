from math import *
def plane_equastion(P, Q, R):
    kx, ky, kz = map(float, [(Q[i] - P[i]) for i in'xyz']) #координаты вектора PQ
    mx, my, mz = map(float, [(R[i] - P[i]) for i in'xyz']) #координаты вектора PR
    px, py, pz = map(float, [P[i] for i in'xyz']) # координаты точки P
    # коэффициенты плоскости заданной в виде ax + by + cz + d = 0
    a = ky*mz - my*kz
    b = mx*kz - kx*mz
    c = kx*my - mx*ky
    d = (-1) * (a*px + b*py + c*pz)
    return {'a': a, 'b': b, 'c': c, 'd': d}

def intersection_of_plain_and_ray(plain, ray_start, ray_vec):
    mx, my, mz = map(float, [ray_start[i] for i in 'xyz']) # координаты начала луча
    vx, vy, vz = map(float, [ray_vec[i] for i in 'xyz']) # координаты направляющего вектора луча
    a, b, c, d = map(float, [plain[i]for i in 'abcd']) # коэффициенты, задающие плоскость
    chislitel = ((a*mx + b*my + c*mz + d) * (-1))
    znamenatel = (a*vx + b*vy + c*vz)
    if znamenatel == 0 or (chislitel/znamenatel) < 0:
        return 0
    alpha = chislitel / znamenatel # коэффициент направляющего вектора луча
    return {'x': (mx + alpha*vx), 'y': (my + alpha*vy), 'z': (mz + alpha*vz)}

def distance_between_dots(A, B):
    return sum([(A[i] - B[i])**2 for i in 'xyz'])**0.5

def rectangle_contains_dot(A, B, C, P):
    # определение точки при прямом угле
    AB = distance_between_dots(A, B)
    AC = distance_between_dots(A, C)
    BC = distance_between_dots(B, C)
    if BC == max([AB, BC, AC]):
        dot_straight_angle, dot_A, dot_B = A, B, C
    elif AC == max([AB, BC, AC]):
        dot_straight_angle, dot_A, dot_B = B, A, C
    else:
        dot_straight_angle, dot_A, dot_B = C, A, B
    # раскладываем радиус-вектор точки по базису векторов AB и AC
    ab_x, ab_y, ab_z = map(float, [(dot_A[i] - dot_straight_angle[i]) for i in 'xyz'])
    ac_x, ac_y, ac_z = map(float, [(dot_B[i] - dot_straight_angle[i]) for i in 'xyz'])
    p_x, p_y, p_z = map(float, [(P[i] - dot_straight_angle[i]) for i in 'xyz'])
    betta_chisliteli = [(p_y*ab_x - p_x*ab_y), (p_y*ab_z - p_z*ab_y), (p_z*ab_x - p_x*ab_z)]
    betta_znamenateli = [(ab_x*ac_y - ab_y*ac_x), (ab_z*ac_y - ab_y*ac_z), (ab_z*ac_x - ab_x*ac_y)]
    for i in range(3):
        if betta_znamenateli[i] != 0:
            betta = betta_chisliteli[i] / betta_znamenateli[i]
            #break
    alpha_chisliteli = [(p_x - betta * ac_x), (p_y - betta * ac_y), (p_z - betta * ac_z)]
    alpha_znamenateli = [ab_x, ab_y, ab_z]
    for i in range(3):
        if alpha_znamenateli[i] != 0:
            alpha = alpha_chisliteli[i] / alpha_znamenateli[i]
            break
    # проверка принадлежности прямоугольнику
    return (0 <= alpha <= 1) and (0 <= betta <= 1)

def reflected_vector(P, Q, R, ray):
    mirror_plain = plane_equastion(P, Q, R)
    n_x, n_y, n_z = map(float, [mirror_plain[i] for i in 'abc'])
    n_len = (n_x**2 + n_y**2 + n_z**2)**0.5
    a_x, a_y, a_z = map(float, [ray[i] for i in 'xyz'])
    a_len = (a_x**2 + a_y**2 + a_z**2)**0.5
    cos_a = (n_x*a_x + n_y*a_y + n_z*a_z) / (n_len*a_len)
    if cos_a > 0:
        n_x, n_y, n_z = -n_x, -n_y, -n_z
    cos_a = (n_x * a_x + n_y * a_y + n_z * a_z) / (n_len * a_len)
    cos_b = 1 - 2*(cos_a**2)
    x = ((2*(a_len**2) - 2*(a_len**2)*cos_b) / (n_len**2)) ** 0.5
    return {'x': (n_x * x + a_x), 'y': (n_y * x + a_y), 'z': (n_z * x + a_z)}

f =  open('input.txt', 'r')
A = dict(zip(['x', 'y', 'z'], [float(x) for x in f.readline().split()]))
B = dict(zip(['x', 'y', 'z'], [float(x) for x in f.readline().split()]))
C = dict(zip(['x', 'y', 'z'], [float(x) for x in f.readline().split()]))
D = dict(zip(['x', 'y', 'z'], [float(x) for x in f.readline().split()]))
P1 = {'x': B['x'], 'y': B['y'], 'z': D['z']}
P2 = {'x': A['x'], 'y': A['y'], 'z': D['z']}
P3 = {'x': A['x'], 'y': D['y'], 'z': D['z']}
P4 = {'x': A['x'], 'y': D['y'], 'z': A['z']}
walls = [
    plane_equastion(A, B, C),
    plane_equastion(P1, P2, D),
    plane_equastion(A, B, P1),
    plane_equastion(A, P2, P3),
    plane_equastion(C, D, P3),
    plane_equastion(B, C, D)
]
V = dict(zip(['x', 'y', 'z'], [float(x) for x in f.readline().split()]))
M = dict(zip(['x', 'y', 'z'], [float(x) for x in f.readline().split()]))
e0 = int(f.readline())
k = int(f.readline())
mirrors = []
for i in range(k):
    P = dict(zip(['x', 'y', 'z'], [float(x) for x in f.readline().split()]))
    Q = dict(zip(['x', 'y', 'z'], [float(x) for x in f.readline().split()]))
    R = dict(zip(['x', 'y', 'z'], [float(x) for x in f.readline().split()]))
    mirrors.append([P, Q, R])

f = open('output.txt', 'w')
k = e0
for i in range(k):
    # точка пересечения со стенками куба
    nearest_wall = {'dist': 10**27, 'coord': 0}
    for wall in walls:
        intersection_coord = intersection_of_plain_and_ray(wall, M, V)
        if intersection_coord != 0:
            if distance_between_dots(intersection_coord, M) < nearest_wall['dist'] and distance_between_dots(intersection_coord, M) != 0:
                nearest_wall['dist'] = distance_between_dots(intersection_coord, M)
                nearest_wall['coord'] = intersection_coord
    intersection_with_mirror_flag = False
    nearest_mirror = {'dist': 10 ** 27, 'coord': 0, 'dots': []}
    for mirror in mirrors:
        P = mirror[0]
        Q = mirror[1]
        R = mirror[2]
        mirror_plain = plane_equastion(P, Q, R)
        mirror_intersection = intersection_of_plain_and_ray(mirror_plain, M, V)
        if mirror_intersection != 0: #луч пересекает плоскость зеркала
            if rectangle_contains_dot(P, Q, R, mirror_intersection):#луч попадает в прямоугльник зеркала
                intersection_with_mirror_flag = True
                mirror_distance = distance_between_dots(M, mirror_intersection)
                if mirror_distance < nearest_mirror['dist'] and mirror_distance != 0:
                    nearest_mirror['dist'] = mirror_distance
                    nearest_mirror['coord'] = mirror_intersection
                    nearest_mirror['dots'] = [P, Q, R]
    if nearest_wall['dist'] < nearest_mirror['dist'] or not intersection_with_mirror_flag:
        f.write('1' + '\n')
        f.write(str(e0) + '\n')
        f.write(str(round(nearest_wall['coord']['x'], 2))+' '+str(round(nearest_wall['coord']['y'], 2))+' '+str(round(nearest_wall['coord']['z'], 2)) + '\n')
        f.write(str(round(V['x'], 2))+' '+str(round(V['y'], 2))+' '+str(round(V['z'], 2)))
        break
    if e0 == 1:
        f.write('0' + '\n')
        f.write(str(nearest_mirror['coord']['x'])+' '+str(nearest_mirror['coord']['y'])+' '+str(nearest_mirror['coord']['z']))
        break
    M = nearest_mirror['coord']
    P = nearest_mirror['dots'][0]
    Q = nearest_mirror['dots'][1]
    R = nearest_mirror['dots'][2]
    V = reflected_vector(P, Q, R, V)
    e0 -= 1




