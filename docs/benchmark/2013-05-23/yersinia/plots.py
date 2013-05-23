from pylab import *

yersinia_read_len = [100, 200, 500, 1000, 2000, 3000]

yersinia01_time_bwamem   = [0.941, 2.090, 10.188, 26.870, 65.722, 100.293]
yersinia01_time_bwasw    = [4.252, 5.382, 13.512, 26.703, 55.096, 90.979]
yersinia01_time_seqalto  = [1.660, 2.117, 10.868, 63.031]
yersinia01_time_snap     = [0.739, 0.789, 0.627, 0.587, 0.739]
yersinia01_time_lisa_lis = [1.681, 2.804, 6.263, 12.178, 24.338, 36.376]
yersinia01_time_lisa_cov = [2.280, 4.272, 10.615, 23.105, 50.501, 88.408]

p1,=plot(yersinia_read_len, yersinia01_time_bwamem)
p2,=plot(yersinia_read_len, yersinia01_time_bwasw)
p3,=plot(yersinia_read_len[0:4], yersinia01_time_seqalto)
p4,=plot(yersinia_read_len[0:5], yersinia01_time_snap)
p5,=plot(yersinia_read_len, yersinia01_time_lisa_lis)
p6,=plot(yersinia_read_len, yersinia01_time_lisa_cov)

legend([p1,p2,p3,p4,p5,p6], 
       ['BWA-MEM', 'BWA-SW', 'SeqAlto', 'SNAP', 'Lisa-LIS', 'Lisa-COV'], loc=2)
xlabel('Duljina ocitanja')
ylabel('Vrijeme izvrsavanja [s]')
show()


yersinia05_time_bwamem   = [1.354, 4.013, 17.605, 40.603, 76.394, 110.365]
yersinia05_time_bwasw    = [4.016, 12.730, 26.367, 56.449, 91.092]
yersinia05_time_seqalto  = [2.194, 5.586, 66.372, 263.149]
yersinia05_time_snap     = [1.135, 1.341, 1.663]
yersinia05_time_lisa_lis = [1.631, 2.711, 5.960, 11.456, 22.565, 33.547]
yersinia05_time_lisa_cov = [1.867, 3.320, 7.984, 17.341, 38.199, 68.602]

p1,=plot(yersinia_read_len, yersinia05_time_bwamem)
p2,=plot(yersinia_read_len[0:5], yersinia05_time_bwasw)
p3,=plot(yersinia_read_len[0:4], yersinia05_time_seqalto)
p4,=plot(yersinia_read_len[0:3], yersinia05_time_snap)
p5,=plot(yersinia_read_len, yersinia05_time_lisa_lis)
p6,=plot(yersinia_read_len, yersinia05_time_lisa_cov)

legend([p1,p2,p3,p4,p5,p6], 
       ['BWA-MEM', 'BWA-SW', 'SeqAlto', 'SNAP', 'Lisa-LIS', 'Lisa-COV'], loc=2)
xlabel('Duljina ocitanja')
ylabel('Vrijeme izvrsavanja [s]')
show()

yersinia10_time_bwamem   = [1.227, 3.745, 15.443, 35.257, 68.372, 99.946]
yersinia10_time_bwasw    = [3.845, 6.092, 14.004, 27.404, 56.873, 88.603]
yersinia10_time_seqalto  = [5.597, 20.619, 86.824, 49.785, 183.418]
yersinia10_time_snap     = [0.930]
yersinia10_time_lisa_lis = [1.510, 2.621, 5.784, 11.063, 21.597, 32.086]
yersinia10_time_lisa_cov = [1.597, 2.857, 6.853, 15.125, 31.159, 57.321]

p1,=plot(yersinia_read_len, yersinia10_time_bwamem)
p2,=plot(yersinia_read_len, yersinia10_time_bwasw)
p3,=plot(yersinia_read_len[0:5], yersinia10_time_seqalto)
p4,=plot(yersinia_read_len[0:1], yersinia10_time_snap)
p5,=plot(yersinia_read_len, yersinia10_time_lisa_lis)
p6,=plot(yersinia_read_len, yersinia10_time_lisa_cov)

legend([p1,p2,p3,p4,p5,p6], 
       ['BWA-MEM', 'BWA-SW', 'SeqAlto', 'SNAP', 'Lisa-LIS', 'Lisa-COV'], loc=2)
xlabel('Duljina ocitanja')
ylabel('Vrijeme izvrsavanja [s]')
show()


