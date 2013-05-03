from pylab import *

# E. Coli
ecoli_read_len = [100, 200, 500]

ecoli_time_bwamem  = [1.059, 2.590, 10.826]
ecoli_time_bwasw   = [3.730, 5.021, 12.355]
ecoli_time_seqalto = [1.623, 1.957, 9.683]
ecoli_time_snap    = [0.689, 0.767, 0.579]
ecoli_time_lisa    = [1.711, 2.928, 6.723]

p1,=plot(ecoli_read_len, ecoli_time_bwamem)
p2,=plot(ecoli_read_len, ecoli_time_bwasw)
p3,=plot(ecoli_read_len, ecoli_time_seqalto)
p4,=plot(ecoli_read_len, ecoli_time_snap)
p5,=plot(ecoli_read_len, ecoli_time_lisa)

legend([p1,p2,p3,p4,p5], ['BWA-MEM', 'BWA-SW', 'SeqAlto', 'SNAP', 'LISA'], loc=2)
xlabel('Duljina ocitanja')
ylabel('Vrijeme izvrsavanja [s]')
show()


# Salmonella
salmonella_read_len = [100, 200, 500, 1000, 2000]
salmonella_time_bwamem  = [1.294, 2.710, 10.430, 37.919, 66.630]
salmonella_time_bwasw   = [4.014, 5.205, 11.751, 24.117, 54.542]
salmonella_time_seqalto = [1.406, 1.878, 8.084, 50.932, 372.40]
salmonella_time_snap    = [0.696, 0.728, 0.757, 1.075]
salmonella_time_lisa    = [1.690, 2.982, 6.712, 12.991, 25.150]

p1,=plot(salmonella_read_len, salmonella_time_bwamem)
p2,=plot(salmonella_read_len, salmonella_time_bwasw)
p3,=plot(salmonella_read_len[0:4], salmonella_time_seqalto[0:4])
p4,=plot([100, 200, 500, 1000], salmonella_time_snap)
p5,=plot(salmonella_read_len, salmonella_time_lisa)

#legend([p1,p2,p4,p5], ['BWA-MEM', 'BWA-SW', 'SNAP', 'LISA'], loc=2)
legend([p1,p2,p3,p4,p5], ['BWA-MEM', 'BWA-SW', 'SeqAlto', 'SNAP', 'LISA'], loc=2)
xlabel('Duljina ocitanja')
ylabel('Vrijeme izvrsavanja [s]')
show()


# Streptococcus
strept_read_len = [100, 200, 500]

strept_time_bwamem  = [1.155, 2.578, 10.691];
strept_time_bwasw   = [3.912, 5.449, 11.662];
strept_time_seqalto = [1.653, 1.747, 8.463];
strept_time_snap    = [0.438, 0.658, 0.490];
strept_time_lisa    = [1.397, 2.342, 4.911];

p1,=plot(strept_read_len, strept_time_bwamem)
p2,=plot(strept_read_len, strept_time_bwasw)
p3,=plot(strept_read_len, strept_time_seqalto)
p4,=plot(strept_read_len, strept_time_snap)
p5,=plot(strept_read_len, strept_time_lisa)


legend([p1,p2,p3,p4,p5], ['BWA-MEM', 'BWA-SW', 'SeqAlto', 'SNAP', 'LISA'], loc=2)
xlabel('Duljina ocitanja')
ylabel('Vrijeme izvrsavanja [s]')
show()



# Yersinia pestis
yersinia_read_len = [100, 200, 500, 1000, 2000, 5000]
yersinia_time_bwamem  = [1.073, 3.841, 16.941, 57.090, 90.931, 191.490]
yersinia_time_bwasw   = [3.855, 5.728, 12.941, 27.229, 59.720, 166.228]
yersinia_time_seqalto = [0.931, 2.145, 11.535, 89.092, 429, 1370]
yersinia_time_snap    = [0.705, 0.740, 1.005, 1.146]
yersinia_time_lisa    = [1.734, 2.989, 6.753, 13.062, 25.648, 66.738]

p1,=plot(yersinia_read_len, yersinia_time_bwamem)
p2,=plot(yersinia_read_len, yersinia_time_bwasw)
p3,=plot(yersinia_read_len[0:4], yersinia_time_seqalto[0:4])
p4,=plot([100, 200, 500, 1000], yersinia_time_snap)
p5,=plot(yersinia_read_len, yersinia_time_lisa)

legend([p1,p2,p3,p4,p5], ['BWA-MEM', 'BWA-SW', 'SeqAlto', 'SNAP', 'LISA'], loc=2)
xlabel('Duljina ocitanja')
ylabel('Vrijeme izvrsavanja [s]')
show()
