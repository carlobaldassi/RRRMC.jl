set logscale
set style data l
N = 1001
tau(t) = 10**t
m(k) = 2**k
set xlabel "MCS"
set ylabel "energy / N"
p for [t=3:3] for [k=5:12] sprintf("test_QSA_rrr_tau%i_seed7001.ALT.M%i.S.txt", m(k)*tau(t), m(k)) u ($1/(N*m(k))):($3) t sprintf("%i", m(k))
