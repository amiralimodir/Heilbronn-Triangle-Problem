import random
import matplotlib.pyplot as plt

def mindist(xy,x2,y2):
    mindist=2**0.5
    for (x1,y1) in xy:
        dist=((x1-x2)**2+(y1-y2)**2)**0.5
        if dist <= mindist:
            mindist = dist
    return mindist


def plot_solution(xy):
    optimal_x=[]
    optimal_y=[]
    for (x,y) in xy:
        optimal_x.append(x)
        optimal_y.append(y)
    plt.figure(figsize=(10, 10))
    plt.scatter(optimal_x, optimal_y, c='red')

    n = len(optimal_x)
    for i in range(n):
        plt.annotate(f"{i}", (optimal_x[i], optimal_y[i]), textcoords="offset points", xytext=(5, 5), ha='center')
        
    x=optimal_x
    y=optimal_y
    minarea=100
    
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                area =abs( 0.5 * (x[i] * (y[j] - y[k]) + x[j] * (y[k] - y[i]) + x[k] * (y[i] - y[j])) )
                if area<=minarea:
                    minarea=area

    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                area =abs( 0.5 * (x[i] * (y[j] - y[k]) + x[j] * (y[k] - y[i]) + x[k] * (y[i] - y[j])) )
                if(area == minarea):
                    trianglex=[x[i],x[j],x[k]]
                    triangley=[y[i],y[j],y[k]]
                    for t in range(3):
                        plt.plot(trianglex, triangley, 'g-')
                        
                    plt.fill(trianglex, triangley)
                    print(area)
                plt.plot([optimal_x[i], optimal_x[j]], [optimal_y[i], optimal_y[j]], 'b-')
                plt.plot([optimal_x[j], optimal_x[k]], [optimal_y[j], optimal_y[k]], 'b-')
                plt.plot([optimal_x[k], optimal_x[i]], [optimal_y[k], optimal_y[i]], 'b-')

    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Optimal Points and Triangles')
    plt.grid(True)
    #plt.savefig('result.jpg')
    plt.show()



ans=[]
ansb=[]
minb=100
maxb=0

for i in range(1000000):
    xy=[]
    xy.append((0,random.uniform(0,1)))
    xy.append((1,random.uniform(0,1)))

    x2=random.uniform(0,1)
    while mindist(xy,x2,0) < 0.167718:
        x2=random.uniform(0,1)
    xy.append((x2,0))

    x2=random.uniform(0,1)
    while mindist(xy,x2,0) < 0.167718:
        x2=random.uniform(0,1)
    xy.append((x2,0))

    x2=random.uniform(0,1)
    while mindist(xy,x2,1) < 0.167718:
        x2=random.uniform(0,1)
    xy.append((x2,1))

    for i in range(2):
        x2=random.uniform(0,1)
        y2=random.uniform(0,1)
        while mindist(xy,x2,y2) < 0.167718:
            x2=random.uniform(0,1)
            y2=random.uniform(0,1)
        xy.append((x2,y2))
    
    minarea=1
    sumb=0
    n=len(xy)
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                area = 0.5 * (xy[i][0] * (xy[j][1] - xy[k][1]) + xy[j][0] * (xy[k][1] - xy[i][1]) + xy[k][0] * (xy[i][1] - xy[j][1])) 
                if area < 0 :
                    b=0
                    area = -1 * area 
                else:
                    b=1
                if area<=minarea:
                    minarea=area
                sumb+=b
    if sumb>=maxb:
        maxb = sumb
        xymax = xy
    if sumb<= minb:
        minb = sumb
        xymin = xy
    ans.append(minarea)
    ansb.append(sumb)

print('max b: ',maxb)
plot_solution(xymax)

print('min b: ', minb)
plot_solution(xymin)

print(sum(ansb)/len(ansb))
print(max(ans))