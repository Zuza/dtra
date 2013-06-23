from pylab import *

yersinia_read_len = [100, 200, 500, 1000, 2000, 3000]
yersinia_read_len_lisa = [100, 200, 500, 1000, 1500, 2000, 3000]

yersinia01_time_bwamem   = [0.941, 2.090, 10.188, 26.870, 65.722, 100.293]
yersinia01_time_bwasw    = [4.252, 5.382, 13.512, 26.703, 55.096, 90.979]
yersinia01_time_seqalto  = [1.660, 2.117, 10.868, 63.031]
yersinia01_time_snap     = [0.739, 0.789, 0.627, 0.587, 0.739]
yersinia01_time_lisa_lis = [3.14, 4.26, 7.36, 12.86, 21.16, 26.97, 49.38]
yersinia01_time_lisa_cov = [2.96, 4.35, 8.76, 19.17, 39.62, 50.61, 84.11]

figure()
p1,=plot(yersinia_read_len, yersinia01_time_bwamem)
p2,=plot(yersinia_read_len, yersinia01_time_bwasw)
p3,=plot(yersinia_read_len[0:4], yersinia01_time_seqalto)
p4,=plot(yersinia_read_len[0:5], yersinia01_time_snap)
p5,=plot(yersinia_read_len_lisa, yersinia01_time_lisa_lis)
p6,=plot(yersinia_read_len_lisa, yersinia01_time_lisa_cov)

legend([p1,p2,p3,p4,p5,p6], 
       ['BWA-MEM', 'BWA-SW', 'SeqAlto', 'SNAP', 'Lisa-LIS', 'Lisa-COV'], loc=2)
xlabel('Read length')
ylabel('Time [s]')
show()


yersinia05_time_bwamem   = [1.354, 4.013, 17.605, 40.603, 76.394, 110.365]
yersinia05_time_bwasw    = [4.016, 12.730, 26.367, 56.449, 91.092]
yersinia05_time_seqalto  = [2.194, 5.586, 66.372, 263.149]
yersinia05_time_snap     = [1.135, 1.341, 1.663]
yersinia05_time_lisa_lis = [2.94, 3.97, 6.69, 12.69, 20.80, 26.94, 47.58]
yersinia05_time_lisa_cov = [3.13, 4.45, 6.90, 13.62, 23.01, 33.86, 78.05]

figure()
p1,=plot(yersinia_read_len, yersinia05_time_bwamem)
p2,=plot(yersinia_read_len[0:5], yersinia05_time_bwasw)
p3,=plot(yersinia_read_len[0:4], yersinia05_time_seqalto)
p4,=plot(yersinia_read_len[0:3], yersinia05_time_snap)
p5,=plot(yersinia_read_len_lisa, yersinia05_time_lisa_lis)
p6,=plot(yersinia_read_len_lisa, yersinia05_time_lisa_cov)

legend([p1,p2,p3,p4,p5,p6], 
       ['BWA-MEM', 'BWA-SW', 'SeqAlto', 'SNAP', 'Lisa-LIS', 'Lisa-COV'], loc=2)
xlabel('Read length')
ylabel('Time [s]')
show()

yersinia10_time_bwamem   = [1.227, 3.745, 15.443, 35.257, 68.372, 99.946]
yersinia10_time_bwasw    = [3.845, 6.092, 14.004, 27.404, 56.873, 88.603]
yersinia10_time_seqalto  = [5.597, 20.619, 86.824, 49.785, 183.418]
#yersinia10_time_snap     = [0.930]
yersinia10_time_lisa_lis = [2.43, 3.63, 6.85, 13.41, 21.98, 28.58, 61.01]
yersinia10_time_lisa_cov = [2.47, 3.65, 6.84, 13.38, 22.09, 28.75, 61.86]

figure()
p1,=plot(yersinia_read_len, yersinia10_time_bwamem)
p2,=plot(yersinia_read_len, yersinia10_time_bwasw)
p3,=plot(yersinia_read_len[0:5], yersinia10_time_seqalto)
#p4,=plot(yersinia_read_len[0:1], yersinia10_time_snap)
p4,=plot(yersinia_read_len_lisa, yersinia10_time_lisa_lis)
p5,=plot(yersinia_read_len_lisa, yersinia10_time_lisa_cov)

legend([p1,p2,p3,p4,p5], 
       ['BWA-MEM', 'BWA-SW', 'SeqAlto', 'Lisa-LIS', 'Lisa-COV'], loc=2)
xlabel('Read length')
ylabel('Time [s]')
show()


