import numpy as np #数値計算ライブラリ
import matplotlib.pyplot as plt #グラフ描画ライブラリ
from scipy.integrate import ode #常微分方程式ライブラリ
import matplotlib.patches as patches


B0 = -0.5 #磁場の強さ[T]
q = 1.0 #電気素量
m= 0.51 #eの質量[Mev/c^2]
R=0.0125

#関数の定義 運動方程式の定義
def newton(t, Y, q, m, B0):
    """
    運動方程式に基づいて状態ベクトルYの微分を計算します。
    Yは状態ベクトルであり、(x, y, z, u, v, w) === (位置, 速度)を含む。
    戻り値はdY/dtです。
    """
    x, y, z = Y[0], Y[1], Y[2]
    u, v, w = Y[3], Y[4], Y[5]
    if x**2+y**2<R**2:
        alpha = 300*q*B0
    else :
        alpha = 0
    

    return np.array([u, v, w,  alpha * v,  -alpha * u ,0]) ## (z,x,vy)になってます
    # 関数の戻り値は、x,y,z,u,v,wの6つの値を持つベクトル

# ODEソルバーの初期設定
r = ode(newton).set_integrator('dopri5')

focus_positions = []

# ラベルをすでに追加したかどうかを追跡するフラグを初期化します。
labels_added = {1: False, 2: False, 3: False, 4: False}

for j in range(1,4): ##繰り返しの回数は3回 運動量を変えていく
    for i in [-0.02,-0.01,0,0.01,0.02]: ## 初期角度を変えていく

        p = 0.3 * j + 1.1  #初期運動量[Mev/c]
        theta = i * np.pi + np.pi/2 #初期角度
        px =  p * np.cos(theta) #初期運動量x成分
        py = p * np.sin(theta) #初期運動量y成分

        rho = p / (300 * B0)  # P[GeV/c]=0.3*B[T]*r[m]から曲率半径の導出

        # 軌跡の座標を計算
        y1 = -2 * R * (R**2 - rho**2) / (3 * R**2 - rho**2)
        x1 = 4 * rho * R**2 / (3 * R**2 - rho**2)
        focus_positions.append([x1,y1])


        # 初期条件の設定
        t0 = 0
        x0 = np.array([0, -2*R , 0])
        v0 = np.array([px, py, 0]) ## (z,x,y)になってます
        initial_conditions = np.concatenate((x0, v0))

        r.set_initial_value(initial_conditions, t0).set_f_params(q, m, B0)


        # 位置を格納するリスト
        positions = []
        positions.append(initial_conditions[:3])
        if j == 1:    
            t1 = 0.05
        if j == 2:
            t1 = 0.04
        if j == 3:
            t1 = 0.035
        
        
        
        dt = 0.001
        while r.successful() and r.t < t1:
            r.integrate(r.t+dt)
            positions.append(r.y[:3]) # 速度ではなく位置だけを保持

        positions = np.array(positions)

        # 結果をプロット
        if j == 1:
            if not labels_added[1]:
                plt.plot(positions[:, 0], positions[:, 1], color='green', label=f'{np.round(p,3)}Mev/c')
                labels_added[1] = True
            else:
                plt.plot(positions[:, 0], positions[:, 1], color='green')
        elif j == 2:
            if not labels_added[2]:
                plt.plot(positions[:, 0], positions[:, 1], color='yellow', label=f'{np.round(p,3)}Mev/c')
                labels_added[2] = True
            else:
                plt.plot(positions[:, 0], positions[:, 1], color='yellow')
        elif j == 3:
            if not labels_added[3]:
                plt.plot(positions[:, 0], positions[:, 1], color='red', label=f'{np.round(p,3)}Mev/c')
                labels_added[3] = True
            else:
                plt.plot(positions[:, 0], positions[:, 1], color='red')

        

        
    

       



## 磁場が存在する範囲を表示
circle = plt.Circle((0, 0), R,fill=True, edgecolor='black',facecolor='lightblue')
plt.gcf().gca().add_artist(circle)
#rectangle=patches.Rectangle((zx,zy),0.5,1,fill=True,edgecolor='black',facecolor='lightblue')
#plt.gcf().gca().add_patch(rectangle)

#焦点の表示
x_focus, y_focus = zip(*focus_positions)
plt.plot(x_focus, y_focus, 'ro', label='focus')


plt.title(f"B={B0}T")
plt.xlabel("X[m]")
plt.ylabel("Z[m]")
plt.legend(fontsize=14)
plt.axis([-0.02, 0.15, -0.2, 0.1])
plt.axis('equal')
plt.grid(True)
plt.show()
