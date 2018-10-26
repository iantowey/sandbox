import csv
import matplotlib.pyplot as plt

m=dict()
with open('/home/ian/Desktop/ACM40640-High-Performance-Comp/Assignment1/code/q2/speedup.csv', 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    for row in spamreader:
        m.setdefault(row[0], []).append(row)
        print ', '.join(row)
        
        
header=m['Threads count '][0]        
del m['Threads count ']
for thrds in sorted(map(int,m.keys())):
    x=thrds
    y=float(m[str(thrds)][0][2])
    plt.plot(x,y,'-o')
    plt.axis([0, 5, 0, 4])
leg = plt.legend(loc='2', ncol=1, shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.xlabel('# Threads')
plt.ylabel('Time')
plt.title('Trapezoid Rule : # Threads v Time')
plt.show()
plt.savefig('/home/ian/Desktop/ACM40640-High-Performance-Comp/Assignment1/code/q2/Speedup.png')
