from pylab import *

banana_read_len = [100, 200, 500, 1000, 1500]#, 2000]

banana02_time_bwamem   = [1.78, 3.37, 11.48, 30.12, 48.41]#, 64.94]
banana02_time_bwasw    = [4.46, 7.01, 14.89, 28.21, 40.88]#, 56.28]
banana02_time_lisa_lis = [42.57, 64.01, 164.23, 316.81, 485.13]#, 661.12]
banana02_time_lisa_cov = [127.98, 221.36, 762.28, 1648.20, 2654.14]#, 3738.15]

figure()
p1,=plot(banana_read_len, banana02_time_bwamem)
p2,=plot(banana_read_len, banana02_time_bwasw)
p3,=plot(banana_read_len, banana02_time_lisa_lis)
p4,=plot(banana_read_len, banana02_time_lisa_cov)

legend([p1,p2,p3,p4], 
       ['BWA-MEM', 'BWA-SW', 'Lisa-LIS', 'Lisa-COV'], loc=2)
xlabel('Read length')
ylabel('Time [s]')
ylim(-100, 3000)
show()

banana05_time_bwamem   = [2.42, 4.82, 15.50, 33.59, 51.51]#, 69.01]
banana05_time_bwasw    = [4.13, 7.24, 15.06, 28.30, 40.49]#, 55.33]
banana05_time_lisa_lis = [32.78, 46.56, 105.37, 197.80, 309.30]#, 405.60]
banana05_time_lisa_cov = [66.04, 117.55, 345.99, 727.36, 1298.14]#, 1749.92]

figure()
p1,=plot(banana_read_len, banana05_time_bwamem)
p2,=plot(banana_read_len, banana05_time_bwasw)
p3,=plot(banana_read_len, banana05_time_lisa_lis)
p4,=plot(banana_read_len, banana05_time_lisa_cov)

legend([p1,p2,p3,p4], 
       ['BWA-MEM', 'BWA-SW', 'Lisa-LIS', 'Lisa-COV'], loc=2)
xlabel('Read length')
ylabel('Time [s]')
show()

banana10_time_bwamem   = [2.36, 5.36, 15.27, 31.85, 48.38]#, 64.67]
banana10_time_bwasw    = [3.48, 6.53, 17.71, 27.97, 39.45]#, 54.91]
banana10_time_lisa_lis = [25.07, 33.42, 66.77, 114.32, 173.26]#, 225.29]
banana10_time_lisa_cov = [36.30, 49.79, 132.84, 261.66, 401.51]#, 528.08]

figure()
p1,=plot(banana_read_len, banana10_time_bwamem)
p2,=plot(banana_read_len, banana10_time_bwasw)
p3,=plot(banana_read_len, banana10_time_lisa_lis)
p4,=plot(banana_read_len, banana10_time_lisa_cov)

legend([p1,p2,p3,p4], 
       ['BWA-MEM', 'BWA-SW', 'Lisa-LIS', 'Lisa-COV'], loc=2)
xlabel('Read length')
ylabel('Time [s]')
show()

