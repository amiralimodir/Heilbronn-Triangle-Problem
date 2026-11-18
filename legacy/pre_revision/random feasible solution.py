import matplotlib.pyplot as plt
import random
def plot_solution(optimal_x, optimal_y):
    #plt.figure(figsize=(10, 10))
    #plt.scatter(optimal_x, optimal_y, c='red')

    n = len(optimal_x)
    #for i in range(n):
        #plt.annotate(f"{i}", (optimal_x[i], optimal_y[i]), textcoords="offset points", xytext=(5, 5), ha='center')
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
                    #for t in range(3):
                        #plt.plot(trianglex, triangley, 'g-')
                    #plt.fill(trianglex, triangley)
                    ans=[]
                    ans.append(area)
                    ans.append((i,j,k))
                    ans.append((x[i],y[i]))
                    ans.append((x[j],y[j]))
                    ans.append((x[k],y[k]))
                    return ans

                #plt.plot([optimal_x[i], optimal_x[j]], [optimal_y[i], optimal_y[j]], 'b-')
                #plt.plot([optimal_x[j], optimal_x[k]], [optimal_y[j], optimal_y[k]], 'b-')
                #plt.plot([optimal_x[k], optimal_x[i]], [optimal_y[k], optimal_y[i]], 'b-')


    # plt.xlim(0, 1)
    # plt.ylim(0, 1)
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.title('Optimal Points and Triangles')
    # plt.grid(True)
    # plt.show()

answerx=[]
answery=[]
maxans=0
for i in range(10000000):
    x = [random.uniform(0,1) for _ in range(6)]
    y = [random.uniform(0,1) for _ in range(6)]
    x.append(0)
    x.append(1)
    x.append(0)
    x.append(1)
    y.append(0)
    y.append(1)
    y.append(1)
    y.append(0)
    a=plot_solution(x,y)
    if a[0]>=maxans:
        maxans=a[0]
        answerx=x
        answery=y

print(maxans)
print(answerx)
print(answery)


