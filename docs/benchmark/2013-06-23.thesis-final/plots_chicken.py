from pylab import *

chicken_read_len = [100, 200, 500, 1000, 1500, 2000]

chicken02_time_bwamem   = [2.48, 4.09, 11.75, 30.18, 48.10, 65.10]
chicken02_time_bwasw    = [5.38, 6.97, 13.42, 27.07, 41.09, 56.36]
chicken02_time_lisa_lis = [55.62, 72.01, 119.29, 201.80, 298.23, 383.67]
chicken02_time_lisa_cov = [59.77, 79.32, 133.36, 234.27, 350.52, 460.01]

figure()
p1,=plot(chicken_read_len, chicken02_time_bwamem)
p2,=plot(chicken_read_len, chicken02_time_bwasw)
p3,=plot(chicken_read_len, chicken02_time_lisa_lis)
p4,=plot(chicken_read_len, chicken02_time_lisa_cov)

legend([p1,p2,p3,p4], 
       ['BWA-MEM', 'BWA-SW', 'Lisa-LIS', 'Lisa-COV'], loc=2)
xlabel('Read length')
ylabel('Time [s]')
ylim(-50, 500)
show()

chicken05_time_bwamem   = [2.74, 5.41, 16.04, 34.13, 51.90, 70.02]
chicken05_time_bwasw    = [5.45, 7.64, 13.52, 28.16, 40.63, 56.21]
chicken05_time_lisa_lis = [54.88, 69.12, 112.01, 185.18, 274.77, 349.95]
chicken05_time_lisa_cov = [55.95, 71.50, 117.44, 192.92, 297.79, 385.99]

figure()
p1,=plot(chicken_read_len, chicken05_time_bwamem)
p2,=plot(chicken_read_len, chicken05_time_bwasw)
p3,=plot(chicken_read_len, chicken05_time_lisa_lis)
p4,=plot(chicken_read_len, chicken05_time_lisa_cov)

legend([p1,p2,p3,p4], 
       ['BWA-MEM', 'BWA-SW', 'Lisa-LIS', 'Lisa-COV'], loc=2)
xlabel('Read length')
ylabel('Time [s]')
show()

chicken10_time_bwamem   = [2.78, 5.31, 15.72, 32.62, 49.58, 66.56]
chicken10_time_bwasw    = [4.17, 6.95, 18.21, 27.82, 41.95, 56.27]
chicken10_time_lisa_lis = [54.60, 67.99, 106.07, 172.11, 254.71, 323.28]
chicken10_time_lisa_cov = [57.67, 68.06, 106.84, 170.29, 257.10, 325.34]

figure()
p1,=plot(chicken_read_len, chicken10_time_bwamem)
p2,=plot(chicken_read_len, chicken10_time_bwasw)
p3,=plot(chicken_read_len, chicken10_time_lisa_lis)
p4,=plot(chicken_read_len, chicken10_time_lisa_cov)

legend([p1,p2,p3,p4], 
       ['BWA-MEM', 'BWA-SW', 'Lisa-LIS', 'Lisa-COV'], loc=2)
xlabel('Read length')
ylabel('Time [s]')
show()

