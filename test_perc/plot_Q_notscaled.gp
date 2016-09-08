set style data l
set logscale
N = 1001
p for [t=5:5] for [m in "32 64 128 256 512 1024 2048"] sprintf("test_QSA_rrr_tau%i_seed7001.ALT.M%s.txt", 10**t, m) u ($1/N):($3) t m
