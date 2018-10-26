import csv
import matplotlib.pyplot as plt

m=dict()
with open('/home/ian/Desktop/ACM40640-High-Performance-Comp/Assignment1/speedup.csv', 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    for row in spamreader:
        m.setdefault(row[0], []).append(row)
        print ', '.join(row)
        
header=m['dim'][0]        
del m['dim']
for dim in sorted(map(int,m.keys())):
    x=[row[1] for row in m[str(dim)]]
    y=[row[4] for row in m[str(dim)]]
    plt.plot(x,y,'-o',label="dim=%d"%dim)
    plt.axis([0, 45, 0, 45])
leg = plt.legend(loc='2', ncol=1, shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.xlabel('# Threads')
plt.ylabel('Speedup')
plt.title('Relative Speedup')
plt.show()
plt.plot(range(0,46),range(0,46))
plt.savefig('/home/ian/Desktop/ACM40640-High-Performance-Comp/Assignment1/Relative_Speedup.png')

