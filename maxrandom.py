import random
ans=[]
for i in range(100000):
    x = [random.uniform(0,1) for _ in range(2)]
    y = [random.uniform(0,1) for _ in range(2)]
    x.append(0.2)
    x.append(1)
    x.append(0.2)
    x.append(1)
    x.append(0)
    y.append(0)
    y.append(1)
    y.append(1)
    y.append(0)
    y.append(0.7)
    minarea=1
    n=len(x)
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                area =abs( 0.5 * (x[i] * (y[j] - y[k]) + x[j] * (y[k] - y[i]) + x[k] * (y[i] - y[j])) )
                if area<=minarea:
                    minarea=area
    ans.append(minarea)
print(max(ans))