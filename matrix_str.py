import numpy as np
import math

def tensor(mat_a, mat_b):
    new_matrix=[]
    for i_a in range(len(mat_a)):
        for i_b in range(len(mat_b)):
            line=[]
            for j_a in range(len(mat_a[i_a])):
                for j_b in range(len(mat_b[i_b])):
                    if mat_a[i_a][j_a] == '0' or mat_b[i_b][j_b] == '0':
                        line.append('0')
                    elif mat_a[i_a][j_a] == '1':
                        line.append(mat_b[i_b][j_b])
                    elif mat_b[i_b][j_b] == '1':
                        line.append(mat_a[i_a][j_a])
                    else:
                        line.append(mat_a[i_a][j_a] + str(' * ') + mat_b[i_b][j_b])
            new_matrix.append(line)
    return new_matrix

def mul(mat_a, mat_b):
    assert(len(mat_a[0])==len(mat_b))
    row_count = len(mat_a)
    col_count = len(mat_b[0])
    new_matrix = []
    for i in range(row_count):
        line = []
        for j in range(col_count):
            res = str('')
            for j_a in range(len(mat_a[0])):
                if mat_a[i][j_a]=='0' or mat_b[j_a][j]=='0':
                    continue
                elif mat_a[i][j_a] == '1':
                    res += mat_b[j_a][j]
                elif mat_b[j_a][j] == '1':
                    res += mat_a[i][j_a]
                else:
                    res += mat_a[i][j_a] + str(' * ') + mat_b[j_a][j]
                if j_a != (len(mat_a[0]) - 1):# and not(mat_a[i][j_a+1]=='0' or mat_b[j_a+1][j]=='0'):
                    res += str(' + ')
            if res=='':
                res = '0'
            line.append(res)
        new_matrix.append(line)
    return new_matrix

def add(mat_a, mat_b):
    assert(len(mat_a)==len(mat_b) and len(mat_a[0])==len(mat_b[0]))
    new_matrix = []
    for i in range(len(mat_a)):
        line = []
        for j in range(len(mat_a)):
            if mat_a[i][j] == '0':
                line.append(mat_b[i][j])
            elif mat_b[i][j] == '0':
                line.append(mat_a[i][j])
            else:
                line.append(mat_a[i][j] + str(' + ') + mat_b[i][j])
        new_matrix.append(line)
    return new_matrix

def sub(mat_a, mat_b):
    assert(len(mat_a)==len(mat_b) and len(mat_a[0])==len(mat_b[0]))
    new_matrix = []
    for i in range(len(mat_a)):
        line = []
        for j in range(len(mat_a)):
            if mat_a[i][j] == '0':
                line.append(mat_b[i][j])
            elif mat_b[i][j] == '0':
                line.append(mat_a[i][j])
            else:
                line.append(mat_a[i][j] + str(' - ') + mat_b[i][j])
        new_matrix.append(line)
    return new_matrix


def print_mat(mat):
    for row in mat:
        for num in row:
            print(num)# + "\t", end ='')
        print("")

op_zero = [ ['0', '0', '0', '0', '0'],
            ['0', '0', '0', '0', '0'],
            ['0', '0', '0', '0', '0'],
            ['0', '0', '0', '0', '0'],
            ['0', '0', '0', '0', '0']]

op_id = [['1', '0', '0', '0', '0'],
         ['0', '1', '0', '0', '0'],
         ['0', '0', '1', '0', '0'],
         ['0', '0', '0', '1', '0'],
         ['0', '0', '0', '0', '1']]

op_a = [ ['0', 's(1)', '0', '0', '0'],
         ['0', '0', 's(2)', '0', '0'],
         ['0', '0', '0', 's(n-1)', '0'],
         ['0', '0', '0', '0', 's(n)'],
         ['0', '0', '0', '0', '0']]

op_ad = [['0', '0', '0', '0', '0'],
         ['s(1)', '0', '0', '0', '0'],
         ['0', 's(2)', '0', '0', '0'],
         ['0', '0', 's(n-1)', '0', '0'],
         ['0', '0', '0', 's(n)', '0']]

op_n = [['0', '0', '0', '0', '0'],
         ['0', '1', '0', '0', '0'],
         ['0', '0', '2', '0', '0'],
         ['0', '0', '0', '(n-1)', '0'],
         ['0', '0', '0', '0', 'n']]

op_sz = [['m', '0', '0', '0', '0'],
         ['0', 'm-1', '0', '0', '0'],
         ['0', '0', '0', '0', '0'],
         ['0', '0', '0', '-m+1', '0'],
         ['0', '0', '0', '0', '-m']]

op_sp = [ ['0', 'fp(m-1)', '0', '0', '0'],
         ['0', '0', 'fp(m-2)', '0', '0'],
         ['0', '0', '0', 'fp(-m+1)', '0'],
         ['0', '0', '0', '0', 'fp(-m)'],
         ['0', '0', '0', '0', '0']]

op_sm = [['0', '0', '0', '0', '0'],
         ['fm(m)', '0', '0', '0', '0'],
         ['0', 'fm(m-1)', '0', '0', '0'],
         ['0', '0', 'fm(-m+2)', '0', '0'],
         ['0', '0', '0', 'fm(-m+1)', '0']]

op_p = [ ['p00-00', 'p00-01', 'p00-02', 'p00-03', 'p00-04', 'p00-05', 'p00-06', 'p00-07', 'p00-08', 'p00-09', 'p00-10', 'p00-11', 'p00-12', 'p00-13', 'p00-14', 'p00-15', 'p00-16', 'p00-17', 'p00-18', 'p00-19', 'p00-20', 'p00-21', 'p00-22', 'p00-23', 'p00-24'],
         ['p01-00', 'p01-01', 'p01-02', 'p01-03', 'p01-04', 'p01-05', 'p01-06', 'p01-07', 'p01-08', 'p01-09', 'p01-10', 'p01-11', 'p01-12', 'p01-13', 'p01-14', 'p01-15', 'p01-16', 'p01-17', 'p01-18', 'p01-19', 'p01-20', 'p01-21', 'p01-22', 'p01-23', 'p01-24'],
         ['p02-00', 'p02-01', 'p02-02', 'p02-03', 'p02-04', 'p02-05', 'p02-06', 'p02-07', 'p02-08', 'p02-09', 'p02-10', 'p02-11', 'p02-12', 'p02-13', 'p02-14', 'p02-15', 'p02-16', 'p02-17', 'p02-18', 'p02-19', 'p02-20', 'p02-21', 'p02-22', 'p02-23', 'p02-24'],
         ['p03-00', 'p03-01', 'p03-02', 'p03-03', 'p03-04', 'p03-05', 'p03-06', 'p03-07', 'p03-08', 'p03-09', 'p03-10', 'p03-11', 'p03-12', 'p03-13', 'p03-14', 'p03-15', 'p03-16', 'p03-17', 'p03-18', 'p03-19', 'p03-20', 'p03-21', 'p03-22', 'p03-23', 'p03-24'],
         ['p04-00', 'p04-01', 'p04-02', 'p04-03', 'p04-04', 'p04-05', 'p04-06', 'p04-07', 'p04-08', 'p04-09', 'p04-10', 'p04-11', 'p04-12', 'p04-13', 'p04-14', 'p04-15', 'p04-16', 'p04-17', 'p04-18', 'p04-19', 'p04-20', 'p04-21', 'p04-22', 'p04-23', 'p04-24'],
         ['p05-00', 'p05-01', 'p05-02', 'p05-03', 'p05-04', 'p05-05', 'p05-06', 'p05-07', 'p05-08', 'p05-09', 'p05-10', 'p05-11', 'p05-12', 'p05-13', 'p05-14', 'p05-15', 'p05-16', 'p05-17', 'p05-18', 'p05-19', 'p05-20', 'p05-21', 'p05-22', 'p05-23', 'p05-24'],
         ['p06-00', 'p06-01', 'p06-02', 'p06-03', 'p06-04', 'p06-05', 'p06-06', 'p06-07', 'p06-08', 'p06-09', 'p06-10', 'p06-11', 'p06-12', 'p06-13', 'p06-14', 'p06-15', 'p06-16', 'p06-17', 'p06-18', 'p06-19', 'p06-20', 'p06-21', 'p06-22', 'p06-23', 'p06-24'],
         ['p07-00', 'p07-01', 'p07-02', 'p07-03', 'p07-04', 'p07-05', 'p07-06', 'p07-07', 'p07-08', 'p07-09', 'p07-10', 'p07-11', 'p07-12', 'p07-13', 'p07-14', 'p07-15', 'p07-16', 'p07-17', 'p07-18', 'p07-19', 'p07-20', 'p07-21', 'p07-22', 'p07-23', 'p07-24'],
         ['p08-00', 'p08-01', 'p08-02', 'p08-03', 'p08-04', 'p08-05', 'p08-06', 'p08-07', 'p08-08', 'p08-09', 'p08-10', 'p08-11', 'p08-12', 'p08-13', 'p08-14', 'p08-15', 'p08-16', 'p08-17', 'p08-18', 'p08-19', 'p08-20', 'p08-21', 'p08-22', 'p08-23', 'p08-24'],
         ['p09-00', 'p09-01', 'p09-02', 'p09-03', 'p09-04', 'p09-05', 'p09-06', 'p09-07', 'p09-08', 'p09-09', 'p09-10', 'p09-11', 'p09-12', 'p09-13', 'p09-14', 'p09-15', 'p09-16', 'p09-17', 'p09-18', 'p09-19', 'p09-20', 'p09-21', 'p09-22', 'p09-23', 'p09-24'],
         ['p10-00', 'p10-01', 'p10-02', 'p10-03', 'p10-04', 'p10-05', 'p10-06', 'p10-07', 'p10-08', 'p10-09', 'p10-10', 'p10-11', 'p10-12', 'p10-13', 'p10-14', 'p10-15', 'p10-16', 'p10-17', 'p10-18', 'p10-19', 'p10-20', 'p10-21', 'p10-22', 'p10-23', 'p10-24'],
         ['p11-00', 'p11-01', 'p11-02', 'p11-03', 'p11-04', 'p11-05', 'p11-06', 'p11-07', 'p11-08', 'p11-09', 'p11-10', 'p11-11', 'p11-12', 'p11-13', 'p11-14', 'p11-15', 'p11-16', 'p11-17', 'p11-18', 'p11-19', 'p11-20', 'p11-21', 'p11-22', 'p11-23', 'p11-24'],
         ['p12-00', 'p12-01', 'p12-02', 'p12-03', 'p12-04', 'p12-05', 'p12-06', 'p12-07', 'p12-08', 'p12-09', 'p12-10', 'p12-11', 'p12-12', 'p12-13', 'p12-14', 'p12-15', 'p12-16', 'p12-17', 'p12-18', 'p12-19', 'p12-20', 'p12-21', 'p12-22', 'p12-23', 'p12-24'],
         ['p13-00', 'p13-01', 'p13-02', 'p13-03', 'p13-04', 'p13-05', 'p13-06', 'p13-07', 'p13-08', 'p13-09', 'p13-10', 'p13-11', 'p13-12', 'p13-13', 'p13-14', 'p13-15', 'p13-16', 'p13-17', 'p13-18', 'p13-19', 'p13-20', 'p13-21', 'p13-22', 'p13-23', 'p13-24'],
         ['p14-00', 'p14-01', 'p14-02', 'p14-03', 'p14-04', 'p14-05', 'p14-06', 'p14-07', 'p14-08', 'p14-09', 'p14-10', 'p14-11', 'p14-12', 'p14-13', 'p14-14', 'p14-15', 'p14-16', 'p14-17', 'p14-18', 'p14-19', 'p14-20', 'p14-21', 'p14-22', 'p14-23', 'p14-24'],
         ['p15-00', 'p15-01', 'p15-02', 'p15-03', 'p15-04', 'p15-05', 'p15-06', 'p15-07', 'p15-08', 'p15-09', 'p15-10', 'p15-11', 'p15-12', 'p15-13', 'p15-14', 'p15-15', 'p15-16', 'p15-17', 'p15-18', 'p15-19', 'p15-20', 'p15-21', 'p15-22', 'p15-23', 'p15-24'],
         ['p16-00', 'p16-01', 'p16-02', 'p16-03', 'p16-04', 'p16-05', 'p16-06', 'p16-07', 'p16-08', 'p16-09', 'p16-10', 'p16-11', 'p16-12', 'p16-13', 'p16-14', 'p16-15', 'p16-16', 'p16-17', 'p16-18', 'p16-19', 'p16-20', 'p16-21', 'p16-22', 'p16-23', 'p16-24'],
         ['p17-00', 'p17-01', 'p17-02', 'p17-03', 'p17-04', 'p17-05', 'p17-06', 'p17-07', 'p17-08', 'p17-09', 'p17-10', 'p17-11', 'p17-12', 'p17-13', 'p17-14', 'p17-15', 'p17-16', 'p17-17', 'p17-18', 'p17-19', 'p17-20', 'p17-21', 'p17-22', 'p17-23', 'p17-24'],
         ['p18-00', 'p18-01', 'p18-02', 'p18-03', 'p18-04', 'p18-05', 'p18-06', 'p18-07', 'p18-08', 'p18-09', 'p18-10', 'p18-11', 'p18-12', 'p18-13', 'p18-14', 'p18-15', 'p18-16', 'p18-17', 'p18-18', 'p18-19', 'p18-20', 'p18-21', 'p18-22', 'p18-23', 'p18-24'],
         ['p19-00', 'p19-01', 'p19-02', 'p19-03', 'p19-04', 'p19-05', 'p19-06', 'p19-07', 'p19-08', 'p19-09', 'p19-10', 'p19-11', 'p19-12', 'p19-13', 'p19-14', 'p19-15', 'p19-16', 'p19-17', 'p19-18', 'p19-19', 'p19-20', 'p19-21', 'p19-22', 'p19-23', 'p19-24'],
         ['p20-00', 'p20-01', 'p20-02', 'p20-03', 'p20-04', 'p20-05', 'p20-06', 'p20-07', 'p20-08', 'p20-09', 'p20-10', 'p20-11', 'p20-12', 'p20-13', 'p20-14', 'p20-15', 'p20-16', 'p20-17', 'p20-18', 'p20-19', 'p20-20', 'p20-21', 'p20-22', 'p20-23', 'p20-24'],
         ['p21-00', 'p21-01', 'p21-02', 'p21-03', 'p21-04', 'p21-05', 'p21-06', 'p21-07', 'p21-08', 'p21-09', 'p21-10', 'p21-11', 'p21-12', 'p21-13', 'p21-14', 'p21-15', 'p21-16', 'p21-17', 'p21-18', 'p21-19', 'p21-20', 'p21-21', 'p21-22', 'p21-23', 'p21-24'],
         ['p22-00', 'p22-01', 'p22-02', 'p22-03', 'p22-04', 'p22-05', 'p22-06', 'p22-07', 'p22-08', 'p22-09', 'p22-10', 'p22-11', 'p22-12', 'p22-13', 'p22-14', 'p22-15', 'p22-16', 'p22-17', 'p22-18', 'p22-19', 'p22-20', 'p22-21', 'p22-22', 'p22-23', 'p22-24'],
         ['p23-00', 'p23-01', 'p23-02', 'p23-03', 'p23-04', 'p23-05', 'p23-06', 'p23-07', 'p23-08', 'p23-09', 'p23-10', 'p23-11', 'p23-12', 'p23-13', 'p23-14', 'p23-15', 'p23-16', 'p23-17', 'p23-18', 'p23-19', 'p23-20', 'p23-21', 'p23-22', 'p23-23', 'p23-24'],
         ['p24-00', 'p24-01', 'p24-02', 'p24-03', 'p24-04', 'p24-05', 'p24-06', 'p24-07', 'p24-08', 'p24-09', 'p24-10', 'p24-11', 'p24-12', 'p24-13', 'p24-14', 'p24-15', 'p24-16', 'p24-17', 'p24-18', 'p24-19', 'p24-20', 'p24-21', 'p24-22', 'p24-23', 'p24-24']]

op_q = [ ['q00', 'q01', 'q00', 'q02', 'q03'],
         ['q10', 'q11', 'q10', 'q12', 'q13'],
         ['q20', 'q21', 'q20', 'q22', 'q23'],
         ['q30', 'q31', 'q30', 'q32', 'q33'],
         ['q40', 'q41', 'q40', 'q42', 'q43']]

op_pq = tensor(op_p, op_q)

op_h = add(add(tensor(op_n, op_id), tensor(op_id, op_sz)), add(tensor(op_a, op_sp), tensor(op_ad, op_sm)))

print_mat(sub(mul(op_h, op_p), mul(op_p, op_h)))
