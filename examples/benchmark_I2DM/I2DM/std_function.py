from std_parameter import *
import math
##########################################
#         interaction potential          #
##########################################
def coulumb_potential(q1,q2,r) -> float:    #unit:eV
    return CONSTANT_K*q1*q2*CONSTANT_E/(r*CONSTANT_A)

def screen_length(q1,q2) -> float:       #unit:A
    return (0.8854*CONSTANT_AB)/(q1**0.23+q2**0.23)/CONSTANT_A

def screen_function(x) -> float:
    return 0.1818*math.e**(-3.2*x)+0.5099*math.e**(-0.9423*x)+0.2802*math.e**(-0.4028*x)+0.0281*math.e**(-0.2016*x)

def potential(q1,q2,r) ->float:       #unit:eV
    return coulumb_potential(q1,q2,r)*screen_function(r/screen_length(q1,q2))
##########################################
#       get collision parameters         #
##########################################
def g0(q1,q2,Er,p,r) -> float:
    return 1-(p/r)**2-potential(q1,q2,r)/Er

def f_theta(q1,q2,Er,p,r) -> float:
    return 1/(r**2 * math.sqrt(g0(q1,q2,Er,p,r)))

def f_time(q1,q2,Er,p,r) -> float:
    return 1/math.sqrt(g0(q1,q2,Er,p,r))-1/math.sqrt(1-(p/r)**2)

def R_min(q1,q2,Er,p) -> float:
    Cons = CONSTANT_K*q1*q2*CONSTANT_E/(2*Er*CONSTANT_A)
    R_up = (math.sqrt(p**2+Cons**2)+Cons)*1.1
    x1 = p
    x2 = R_up
    r  = (x1+x2)/2
    value = g0(q1,q2,Er,p,r)
    while abs(value)>1e-14 or value==0:
        if value == 0:
            return r
        elif value > 0:
            x2 = r
        else:
            x1 = r
        r = (x1+x2)/2
        value = g0(q1,q2,Er,p,r)
    return x2
def get_theta(q1,q2,Er,p,Rmin,dx,fac) -> float:
    integral_value = 0
    x = Rmin
    while x < 1e10:
        i = 0
        while i < 10 and x < 1e10:
            x += dx/2
            integral_value += dx*f_theta(q1,q2,Er,p,x)
            x += dx/2
            i += 1
        if dx < 1000:
            dx *= fac
        else:
            dx *= 2
    return math.pi - 2 * p * integral_value
def get_tao(q1,q2,Er,p,Rmin,dx,fac) -> float:
    integral_value = 0
    x = Rmin
    while x < 1e10:
        i = 0
        while i < 10 and x < 1e10:
            x += dx/2
            integral_value += dx*f_time(q1,q2,Er,p,x)
            x += dx/2
            i += 1
        if dx < 1000:
            dx *= fac
        else:
            dx *= 2
    return math.sqrt(Rmin**2-p**2)-integral_value
##########################################
#         electronic stopping            #     
##########################################
def get_x_nl(q1,q2,m1,m2,E,beta) -> float:
    Er = m2*E/(m1+m2)
    reduced_E = screen_length(q1,q2)*Er*CONSTANT_A/(CONSTANT_K*q1*q2*CONSTANT_E)
    return beta*reduced_E**0.075
def get_Se(q1,q2,m1,E,alpha) -> float:
    '''
    unit:eV/A^2
    '''
    Em = 0.5*q1**(4/3)*m1*CONSTANT_MU*CONSTANT_VB**2/CONSTANT_EV
    k_LS = 1.212*q1**(7/6)*q2/((q1**(2/3)+q2**(2/3))**1.5*m1**0.5)
    p = E/Em
    factor = 1/((p/(math.log(p+1/p-2+math.e)))**(1.425/2)+(1/p)**(1.425/2))**(1/1.425)
    #print("Em: ", Em, "\n", "k_LS: ", k_LS, "\n", "factor: ", factor)
    return alpha * k_LS * factor * Em**0.5
def electron_energy_loss(q1,q2,m1,m2,E,p,a0,N,alpha,beta):
    #Nonlocal part
    x_nl = get_x_nl(q1,q2,m1,m2,E,beta)
    x_loc = 1- x_nl
    Se = get_Se(q1,q2,m1,E,alpha)
    pmax = a0/2
    #print("pmax: ", pmax)
    a = 1.45*screen_length(q1,q2)/(0.3*q1**0.4)
    nonlocal_part = N*Se*(x_nl+x_loc*(1+pmax/a)*math.exp(-pmax/a))
    local_part = x_loc*Se*math.exp(-p/a)/(2*math.pi*a**2)
    #print("x_nl: ", x_nl, "\n", "x_loc: ", x_loc, "\n", "Se: ", Se, "\n", "Q_nl: ", nonlocal_part, "\n", "Q_loc: ", local_part)
    return nonlocal_part + local_part
####################################
#   vector and matrices operation  $
####################################
def get_norm_of_vector(vector:list)->float:
    return dis_bet_points(vector,[0,0,0])
def get_direction_of_vector(vector:list):
    value = get_norm_of_vector(vector)
    if value == 0 or vector ==1:
        return vector
    else:
        return [x/value for x in vector]
def get_velocity_value_from_energy(energy:float,mass:int) -> float:
    return math.sqrt(2*energy* CONSTANT_EV/(mass*CONSTANT_MU))*1e-5     #unit:A/fs 
def dis_bet_points(p1,p2) -> float:
    '''Compute the distance between two points'''
    return math.sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2+(p1[2]-p2[2])**2)
def dis_bet_line_point(origin,direction,point) -> list:
    '''Compute the distance line and point in space
    the formula of line can be described as:
    (x,y,z)=(a1,a2,a3)+t(r1,r2,r3)
    '''
    x = origin
    r = direction
    p = point
    A = r[0]**2+r[1]**2+r[2]**2
    B = 2*(r[0]*(x[0]-p[0])+r[1]*(x[1]-p[1])+r[2]*(x[2]-p[2]))
    C = (x[0]-p[0])**2+(x[1]-p[1])**2+(x[2]-p[2])**2
    t = -B/(2*A)
    distance_foot_p = math.sqrt(C-B**2/(4*A))     
    foot = [value+t*v for (value,v) in zip(x,r)]   
    distance_foot_x = dis_bet_points(x,foot)  
    return [t,distance_foot_p,distance_foot_x,foot]
def vector_cross(v1:list,v2:list) -> list:
    return [v1[1]*v2[2]-v1[2]*v2[1],
            v1[2]*v2[0]-v1[0]*v2[2],
            v1[0]*v2[1]-v1[1]*v2[0]]
def multi_matrix(m1:list,m2:list) -> list:
    results = []
    for i in range(3):
        result = []
        for j in range(3):
            value = 0
            for k in range(3):
                value += m1[i][k]*m2[k][j]
            result.append(value)
        results.append(result)
    return results
def multi_matrices(matrices:list) -> list:
    '''multipy all matrices'''
    result = [[1.0,0,0],[0,1,0],[0,0,1]]
    for matrix in matrices:
        result=multi_matrix(result,matrix)
    return result
def vector_multi_matrix(vector:list,matrix:list) -> list:
    result = []
    for i in range(3):
        value = 0
        for j in range(3):
            value += vector[j]*matrix[j][i]
        result.append(value)
    return result
def matrix_dot3(m:list) -> float:
    return (m[0][0]*m[1][1]*m[2][2]+
            m[0][1]*m[1][2]*m[2][0]+
            m[1][0]*m[2][1]*m[0][2]-
            m[0][2]*m[1][1]*m[2][0]-
            m[0][0]*m[1][2]*m[2][1]-
            m[2][2]*m[0][1]*m[1][0])
def matrix_trans(m:list) -> list:
    results = []
    for i in range(3):
        result = []
        for j in range(3):
            result.append(m[j][i])
        results.append(result)
    return results
def residual_matrix_dot(x,y,matrix:list) ->float:
    results = []
    for i in range(3):
        if i != x:
            result = []
            for j in range(3):
                if j!=y:
                    result.append(matrix[i][j])
            results.append(result)
    return results[0][0]*results[1][1]-results[0][1]*results[1][0]
def inverse_matrix(matrix:list) -> list:
    dot = matrix_dot3(matrix)
    if dot ==0:
        raise ValueError
    results = []
    matrix_T = matrix_trans(matrix)
    for i in range(3):
        result = []
        for j in range(3):
            result.append((-1)**(i+j)*residual_matrix_dot(i,j,matrix_T)/dot)
        results.append(result)
    return results
def rotation_matrix(axis_type:str,alpha:float) -> list:
    if axis_type == 'x' or axis_type == 'X':
        return [[1.0,0.0,0.0],
                [0.0,math.cos(alpha),math.sin(alpha)],
                [0.0,-math.sin(alpha),math.cos(alpha)]]
    elif axis_type == 'y' or axis_type == 'Y':
        return [[math.cos(alpha),0.0,-math.sin(alpha)],
                [0.0,1.0,0.0],
                [math.sin(alpha),0.0,math.cos(alpha)]]
    else:
        return [[math.cos(alpha),math.sin(alpha),0.0],
                [-math.sin(alpha),math.cos(alpha),0.0],
                [0.0,0.0,1.0]]
def vector_rotation(vector:list,spiale:list,theta:float) -> list:
    '''Compute the vector after rotating around spiale with a angle of theta.
       vector = (v1,v2,v3)
       spiale = (s1,s2,s3) = s_value(cos(alpha)sin(beta),sin(alpha)sin(beta)),cos(beta))
       beta:[0,pai]
       alpha:(-pai,pai]
       theta:(-pai,pai]
    '''
    s_value = dis_bet_points(spiale,[0,0,0])
    s_radius = math.sqrt(spiale[0]**2+spiale[1]**2)
    if s_value ==0:
        raise ValueError
    if s_radius!=0:
        beta = math.acos(spiale[2]/s_value) 
        value_x = spiale[0]/s_radius
        value_y = spiale[1]/s_radius
        if value_y>=0:
            alpha = math.acos(value_x)
        else:
            alpha = -math.acos(value_x)
        final_rotation_matrix = multi_matrices([rotation_matrix('z',-alpha),
                                                rotation_matrix('y',-beta),
                                                rotation_matrix('z',theta),
                                                rotation_matrix('y',beta),
                                                rotation_matrix('z',alpha)])
    else:
        final_rotation_matrix = rotation_matrix('z',theta*spiale[2]/s_value)
    return vector_multi_matrix(vector,final_rotation_matrix)