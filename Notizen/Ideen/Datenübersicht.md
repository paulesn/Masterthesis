## Vergleich visSpanner und Point Spanner (Full Spanner)

Diese Tabelle vergleicht für mehere Größen von theta was passiert wenn man den 0025.32 Graph als theta basis nimmt oder den vollen vis graph für die gleichen Punkte als Basis nimmt

! nochmal anschauen

|               |              |            |     |             |         |     |             |
| ------------- | ------------ | ---------- | --- | ----------- | ------- | --- | ----------- |
| theta         | full spanner |            |     | vis spanner |         |     |             |
| theta         | avg t        | max t      |     | avg t       | max t   |     |             |
| 12            | 1.02291      | 1.16676    |     | 1.02382     | 1.14081 |     | 0.999111172 |
| 24            | 1.00887      | 1.08618    |     | 1.0089      | 1.07866 |     | 0.999970265 |
| 32            | 1.00543      | 1.13718    |     | 1.00543     | 1.13718 |     | 1           |
| 50            | 1.00356      | 1.08021    |     | 1.00357     | 1.03937 |     | 0.999990036 |
| 64            | 1.00235      | 1.05597    |     | 1.00236     | 1.02999 |     | 0.999990024 |
| 75            | 1.0018       | 1.26097    |     | 1.00179     | 1.02885 |     | 1.000009982 |
| 100           | 1.00119      | 1.02147    |     | 1.00119     | 1.01966 |     | 1           |
| 124           | 1.00115      | 2.20951    |     | 1.00113     | 1.01612 |     | 1.000019977 |
| 128           | 1.00102      | 2.62111    |     | 1.00097     | 1.01602 |     | 1.000049952 |
|               |              |            |     |             |         |     |             |
|               |              |            |     |             |         |     |             |
|               |              |            |     |             |         |     |             |
| full spanner: |              | 1466966601 |     |             |         |     |             |
|               |              | 7849912    |     |             |         |     |             |

## Edges
Hier ist ein Plot für 0025.32. Er zeigt für an der x-axis alle Edges der originialen vis-Graphen sortiert nach ihrem t-Value im Theta Spanner mit 24 cones. Die Blaue linie gibt an was der mean t-Value wäre, wenn man alle Kanten nach dieser dem Spanner hinzufügt. Die rote Linie zeigt für diesen Fall den max t-value


![[running_mean_edges.png]]

![[alltheta.png]]
## 0025.32
Diese Tabelle zeigt wie viel % der $7849912$ Edges im 0025.32 Visibility Graph zum $\theta-$Spanner hinzugefügt werden müssen um ein gewisses target maximum t-value zu erreichen:

| $\theta$\\ target $t$ | $1.1$    | $1.075$ | $1.05$     | $1.025$ | $1.01$  | $1.00$  | $1.001$ |
| --------------------- | -------- | ------- | ---------- | ------- | ------- | ------- | ------- |
| 12                    | 1.356%   | 4.506%  | 14.138%    | 43.205% | 68.662% | 78.919% | 89.091% |
| 24                    | 0.005%   | 0.116%  | 1.231%     | 13.268% | 41.413% | 58.356% | 79.747% |
| 32                    | 0.00028% | 0.009%  | 0.228%     | 5.982%  | 30.011% | 47.781% | 73.964  |
| 50                    | 0        | 0       | 0.0045%    | 1.242%  | 14.680% | 32.550% | 63.959% |
| 64                    | 0        | 0       | 0.00028%   | 0.341%  | 8.759%  | 23.805% | 57.032% |
| 75                    | 0        | 0       | 2.547e-05% | 0.0957% | 4.562%  | 16.574% | 51.895% |
| 100                   | 0        | 0       | 2.547e-05% | 0.0064% | 2.568%  | 12.116% | 44.420  |
| 124                   | 0        | 0       | 2.547e-05% | 0.0012% | 1.094%  | 7.799%  | 38.378% |
| 128                   | 0        | 0       | 2.547e-05% | 0.0019% | 0.911%  | 7.135%  | 37.458% |


Die folgende Tabelle zeigt die gesamtgröße in Edges des Spanners, nachdem die extra Edges hinzugefügt wurden um das target t zu erreichen:

| $\theta$\\ target $t$ | **1.1** | **1.075** | **1.05** | **1.025** | **1.01** | **1.005** | **1.001** |
| --------------------- | ------- | --------- | -------- | --------- | -------- | --------- | --------- |
| **12**                | 444857  | 568486    | 946538   | 2087398   | 3086577  | 3489138   | 3888428   |
| **24**                | 639033  | 643371    | 687156   | 1159574   | 2264265  | 2929285   | 3768895   |
| **32**                | 770097  | 770458    | 779052   | 1004878   | 1947992  | 2645503   | 3673168   |
| **50**                | 1012640 | 1012640   | 1012817  | 1061377   | 1588838  | 2290205   | 3523001   |
| **64**                | 1169282 | 1169282   | 1169293  | 1182666   | 1513066  | 2103612   | 3407760   |
| **75**                | 1275094 | 1275094   | 1275095  | 1278851   | 1454139  | 1925612   | 3311949   |
| **100**               | 1505532 | 1505532   | 1505533  | 1505786   | 1606314  | 1981086   | 3248999   |
| **124**               | 1690590 | 1690590   | 1690592  | 1690684   | 1776506  | 2302834   | 4703254   |
| **128**               | 1718876 | 1718876   | 1718877  | 1718951   | 1754628  | 1998940   | 3189082   |

And here is the table in calc with a heatmap:

![[Pasted image 20250924210447.png]]


Hier zu sehen ist die Prozentuale Anzahl der Edges im Verhältnis zum Orginalen Vis Graph 0025.32
![[Pasted image 20251014152512.png]]
## t vs freq of use in paths
Die folgenden Plots zeigen wie häufig eine Kante in den (500.000) zufälligen Pfaden verwendet wurde vs ihren t-wert.

### $\theta = 24$
![[Pasted image 20250930182207.png]]
Pearson correlation: 0.9999814529909805  (p = 0.0 )

### $\theta = 32$
![[Pasted image 20250930182313.png]]
Pearson correlation: 0.9999900115358752  (p = 0.0 )
### $\theta = 50$
![[Pasted image 20250930182639.png]]
Pearson correlation: 0.9999880933137989  (p = 0.0 )
### $\theta = 64$
![[Pasted image 20250930182822.png]]
Pearson correlation: 0.9999905344074431  (p = 0.0 )
### $\theta = 75$
![[Pasted image 20250930182956.png]]
Pearson correlation: 0.9999851903523117  (p = 0.0 )
### $\theta = 100$
![[Pasted image 20250930183100.png]]
Pearson correlation: 0.9999822714012421  (p = 0.0 )
### $\theta = 124$
![[Pasted image 20250930183201.png]]
Pearson correlation: 0.9999849457860321  (p = 0.0 )
### $\theta = 128$
![[Pasted image 20250930183234.png]]
Pearson correlation: 0.999989524778978  (p = 0.0 )



## Alles auf einmal, mit dem mean t Wert für jeden Count

![[Pasted image 20250930181800.png]]
![[Pasted image 20250930181901.png]]
![[Pasted image 20250930182004.png]]
