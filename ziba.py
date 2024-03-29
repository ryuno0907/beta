import numpy as np #数値計算ライブラリ
import matplotlib.pyplot as plt #グラフ描画ライブラリ
from scipy.integrate import ode #常微分方程式ライブラリ
import matplotlib.patches as patches
import math


#B0 = 0.604 #磁場の強さ[T]
#B1 = 0.01
q = 1.0 #電気素量
m= 0.51 #eの質量[Mev/c^2]
R=0.0125
R1=0.0175
a=-0.00
b=-2*R1
#関数の定義 運動方程式の定義
def newton(t, Y, q, m):
    """
    運動方程式に基づいて状態ベクトルYの微分を計算します。
    Yは状態ベクトルであり、(x, y, z, u, v, w) === (位置, 速度)を含む。
    戻り値はdY/dtです。
    """
    x, y, z = Y[0], Y[1], Y[2]
    u, v, w = Y[3], Y[4], Y[5]
    if (y-3.676/2199/2)**2+x**2<0.015**2:
        r=np.sqrt((y-3.676/2199/2)**2+x**2)
        alpha = -300*q*(-2199*r**2+3.676*r+0.6194)
    else :
        alpha = 0
    

    return np.array([u, v, w,  alpha * v,  -alpha * u ,0]) ## (z,x,vy)になってます
    # 関数の戻り値は、x,y,z,u,v,wの6つの値を持つベクトル

# ODEソルバーの初期設定
r = ode(newton).set_integrator('dopri5')



# ラベルをすでに追加したかどうかを追跡するフラグを初期化します。
labels_added = {1: False, 2: False, 3: False, 4: False}

for j in range(1,5): ## 運動量を変えていく
    for i in [-0.02,-0.01,0,0.01,0.02]: ## 初期角度を変えていく

        p = 0.3 * j + 1.1  #初期運動量[Mev/c]
        theta = i * np.pi + np.pi/2 #初期角度
        if a==0:
            alpha=0
        else:
            alpha = np.pi/2-math.atan(b/a)
        px =  p * np.cos(theta-alpha) #初期運動量x成分
        py = p * np.sin(theta-alpha) #初期運動量y成分


        # 初期条件の設定
        t0 = 0
        x0 = np.array([a, b, 0])
        v0 = np.array([px, py, 0]) ## (z,x,y)になってます
        initial_conditions = np.concatenate((x0, v0))

        r.set_initial_value(initial_conditions, t0).set_f_params(q, m)


        # 位置を格納するリスト
        positions = []
        positions.append(initial_conditions[:3])
        if j == 1:    
            t1 = 0.045
        if j == 2:
            t1 = 0.04
        if j == 3:
            t1 = 0.04
        if j == 4:
            t1 = 0.035
        
        
    

        dt = 0.001
        while r.successful() and r.t < t1:
            r.integrate(r.t+dt)
            positions.append(r.y[:3]) # 速度ではなく位置だけを保持

        positions = np.array(positions)

        # 結果をプロット
        
        if j == 1:
            if not labels_added[1]:
                plt.plot(positions[:, 0], positions[:, 1], color='red', label=f'{np.round(p,3)}Mev/c')
                labels_added[1] = True
            else:
                plt.plot(positions[:, 0], positions[:, 1], color='red')
        elif j == 2:
            if not labels_added[2]:
                plt.plot(positions[:, 0], positions[:, 1], color='blue', label=f'{np.round(p,3)}Mev/c')
                labels_added[2] = True
            else:
                plt.plot(positions[:, 0], positions[:, 1], color='blue')
        
        elif j == 3:
            if not labels_added[3]:
                plt.plot(positions[:, 0], positions[:, 1], color='green', label=f'{np.round(p,3)}Mev/c')
                labels_added[3] = True
            else:
                plt.plot(positions[:, 0], positions[:, 1], color='green')
        elif j == 4:
            if not labels_added[4]:
                plt.plot(positions[:, 0], positions[:, 1], color='orange', label=f'{np.round(p,3)}Mev/c')
                labels_added[4] = True
            else:
                plt.plot(positions[:, 0], positions[:, 1], color='orange')
        #print(positions)
        

        
    

       



## 磁場が存在する範囲を表示

a1=0.03
zx=-a1/2

rectangle=patches.Rectangle((zx,zx),a1,a1,fill=True,edgecolor='black')
plt.gcf().gca().add_patch(rectangle)
circle = plt.Circle((0, 0), R,fill=True, edgecolor='black',facecolor='lightblue')
plt.gcf().gca().add_artist(circle)

# タイトルやラベルなどを付ける場合
#plt.title(f"B={B0}T")
plt.xlabel("X[m]")
plt.ylabel("Z[m]")
# 凡例を大きく表示
plt.legend(fontsize=14)
#plt.xlim(-0.0125,0.05)
#plt.axis([-0.1, 0.05, -0.2, 0.1])
##plt.axis([-0.5,0.5,-100,100])
plt.xlim(-0.05,0.05)
plt.ylim(-0.05,0.05)

plt.axis('equal')


plt.grid(True)

#plt.show()
plt.xlabel("X[m]")
plt.ylabel("Y[m]")
plt.axis('equal')
#plt.xlim(0,0.2)
#plt.ylim(0,0.2)
plt.legend()
plt.grid(True)
plt.show()