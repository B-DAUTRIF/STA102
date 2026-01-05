---
title: "Régression linéaire simple"
output: html_notebook
---

Ce bloc-note reprend le contenu du cours STA102 du Conservatoire National des Arts et Métier (CNAM) sur les modéles linéaires. Son but est de décrire les étapes d'une régression linéaire simple, en rassemblant à la fois le matériel théorique et pratique pour travailler efficacement sur des données.

Ce script est destiné à un usage personnel. J'essaie au maximum qu'il puisse être utilisé tel quel pour effectuer les calculs et la production des différents résultats (matrices, vecteurs paramètres et de variance, testes statistiques) : a condition pour cela d'effectuer les changement nécessaires pour l'import de nouveaux jeux de données.

Un éventuel lecteur intéressé par ce document est invité à laisser suggestions et remarques (par le biais de GitHub).

```{r Chargement des librairies}
# install.packages("olsrr")

library(DescTools)
library(lmtest)
library(olsrr)
library(ggplot2)
library(ggfortify)
```

Chargement du jeu de données :

```{r Chargement du jeu de données}
df <- read.delim("D:/Etudes/2025_2026/STA102/Regression_simple/appartements.txt")
```


### Contexte : 

Soit $Y$ une variable **quantitative** (à expliquer) et $X$ une variable explicative :

$Y= \{ y_{1},...,y_{i}...,y_{n} \} $

$X= \{ x_{1},...,x_{i}...,x_{n} \} $

Et le vecteur d'observations composé des couples :

$(x_{i} \ ; \ y_{i}) \quad ; \quad i\in [1\ ;\ n]$

### Modèle de régression simple : 

$$Y=X\beta+\epsilon$$

```{r Modele}
modele=lm(Prix~Superficie, data = df)
```

Représentation graphique :

L'affichage du nuage de points donne une première impression sur l'existence d'une corrélation entre $X$ et $Y$ :

```{r Affchage du nuage de points}
# Graphique :
plot(df[2:3])
```

Quantification de la force de la relation linéaire entre $X$ et $Y$ :

$$r=\frac{\frac{1}{n}\Sigma^{n}_{i=1}(x_{i}-\bar{x})(y_{i}-\bar{y})}{\sqrt{\frac{1}{n}\Sigma^{n}_{i=1}(x_{i}-\bar{x})^{2} \ \frac{1}{n}\Sigma^{n}_{i=1}(y_{i}-\bar{y})^{2}}}$$
$$=\frac{\Sigma^{n}_{i=1}(x_{i}-\bar{x})(y_{i}-\bar{y})}{\sqrt{\Sigma^{n}_{i=1}(x_{i}-\bar{x})^{2} \ \Sigma^{n}_{i=1}(y_{i}-\bar{y})^{2}}}$$
$$=\frac{\Sigma^{n}_{i=1}x_{i}y_{i} \ \  n\bar{x}\bar{y}} {\sqrt{(\Sigma^{n}_{i=1}x_{i}^{2}-n\bar{x}^{2}) \ (\Sigma^{n}_{i=1}y_{i}^{2}-n\bar{y}^{2})}}=\frac{S_{XY}}{S_{X}S_{Y}}$$

Voir théorème de **König-Huygens**.      

Propriétés : 

- $r\in[-1,1]$

**Attention :**

- $X \perp Y \Rightarrow r=0$

- $r=0 \ \nRightarrow X \perp Y$

### Principe des moindres carrés

On cherche à approximer le nuage de points
$(x_{i}, y_{i})$ par une droite d’équation $\hat{y} = a+bx$ de telle sorte que
$\Sigma_{i=1}^{n}(y_{i}− \hat{y}_{i})$ soit minimale :

$$e_{i}=y_{i}-\hat{y_{i}}$$

Avec :

$$Min(\Sigma_{i=1}^{n}y_{i}-\hat{y_{i}})^{2}$$

### Détermination des coefficients : 

Soit : 

$F(a;b)=(\Sigma_{i=1}^{n}y_{i}-a -bx_{i})^{2}$

$$\frac{\partial \ F(a;b)}{\partial \ a}=-2\Sigma_{i=1}^{n}(y_{i}-a -bx_{i})=0$$
$$\Rightarrow \ \Sigma_{i=1}^{n}e_{i}=0$$

$$\frac{\partial \ F(a;b)}{\partial \ b}=-2\Sigma_{i=1}^{n}x_{i}(y_{i}-a -bx_{i})=0$$
$$\Rightarrow \ \Sigma_{i=1}^{n}x_{i}e_{i}=0$$

Donc : 

- $E[\ e\ ]=0$
- $Cov(x,e)=0$

Estimation de **A** :

$$\Sigma_{i=1}^{n}(y_{i}-a -bx_{i})=\Sigma_{i=1}^{n}y_{i} -na -\Sigma^{n}_{i=1}bx_{i}=0$$
$$\Rightarrow \ \frac{\Sigma^{n}_{i=1}y_{i}}{n} -a -b\frac{\Sigma^{n}_{i=1}x_{i}}{n}= \bar{y}-a-b\bar{x}= 0$$
$$\Rightarrow \ a=\bar{y}-b\bar{x}$$

Estimation de **B** :

$$\Sigma^{n}_{i=1}x_{i}(y_{i}-a -bx_{i})=\Sigma^{n}_{i=1}x_{i}(y_{i}-\bar{y}+b\bar{x} -bx_{i})=\Sigma^{n}_{i=1}x_{i}(y_{i}-\bar{y})+b\bar{x}\Sigma^{n}_{i=1}x_{i}-b\Sigma^{n}_{i=1}x_{i}^{2}= 0$$
$$\Sigma^{n}_{i=1}x_{i}(y_{i}-\bar{y})+bn\bar{x}^{2}-b\Sigma^{n}_{i=1}x_{i}^{2}= 0$$
$$b\Sigma^{n}_{i=1}x_{i}^{2}-bn\bar{x}^{2}= \Sigma^{n}_{i=1}x_{i}(y_{i}-\bar{y})$$
$$b(\Sigma^{n}_{i=1}x_{i}^{2}-n\bar{x}^{2})= \Sigma^{n}_{i=1}x_{i}(y_{i}-\bar{y})$$
$$b= \frac{\Sigma^{n}_{i=1}x_{i}(y_{i}-\bar{y})}{\Sigma^{n}_{i=1}x_{i}^{2}-n\bar{x}^{2}}$$

Avec :

$$
x_{i}=x_{i}+\bar{x}-\bar{x} \ \Rightarrow \ 
\sum_{i=1}^n x_i (y_i-\bar y)
= \sum_{i=1}^n [(x_i-\bar x)+\bar x](y_i-\bar y)
$$

$$
= \sum_{i=1}^n (x_i-\bar x)(y_i-\bar y)+ \bar x \sum_{i=1}^n (y_i-\bar y)
$$

et :

$$
\sum_{i=1}^n (y_i-\bar{y})=0 \ \Rightarrow \ \sum_{i=1}^n x_{i} (y_{i}-\bar{y})
= \sum_{i=1}^n (x_{i}-\bar{x})(y_{i}-\bar{y})
$$

$$
(x_i-\bar{x})^2 = x_i^2 - 2\bar x x_i + \bar{x}^2
\quad \Rightarrow \quad \sum_{i=1}^n (x_i-\bar{x})^2
= \sum_{i=1}^n x_i^2-2\bar{x} \sum_{i=1}^n x_i+n\bar{x}^2
= \sum_{i=1}^n x_i^2-2n\bar{x}^{2}+n\bar{x}^2
$$

$$
\sum_{i=1}^n x_i = n\bar x\quad \Rightarrow \quad \sum_{i=1}^n (x_i-\bar x)^2=\sum_{i=1}^n x_i^2 - n\bar x^2
$$

Donc :

$$
b= \frac{\sum_{i=1}^n x_i(y_i-\bar y)}
       {\sum_{i=1}^n x_i^2 - n\bar x^2}
= \frac{\sum_{i=1}^n (x_i-\bar x)(y_i-\bar y)}
       {\sum_{i=1}^n (x_i-\bar x)^2}
$$

On fais le lien avec le coefficient de corrélation : 

$$
r = \frac{s_{XY}}{s_X s_Y} \ \Longrightarrow \ 
s_{XY} = r \times s_X s_Y \ \Longrightarrow \
\frac{s_{XY}}{s_X^2} = \frac{r \times s_X s_Y}{s_X^2} \ \Longrightarrow \ 
\frac{s_{XY}}{s_X^2} = r \frac{s_Y}{s_X}
$$

$$\Rightarrow \ b=\frac{\Sigma^{n}_{i=1}x_{i}y_{i} \ \  n\bar{x}\bar{y}} {\sum_{i=1}^n x_{i}^{2}-n\bar x^{2}}=\frac{S_{XY}}{S_{X}^{2}}=r\frac{s_Y}{s_X}$$

#### Equation de la droite des moindres carrés :

$$\hat{y}=a+bx=\bar{y}-b\bar{x}+r\frac{S_{Y}}{S_{X}}x=\bar{y}-r\frac{S_{Y}}{S_{X}}\bar{x}+r\frac{S_{Y}}{S_{X}}x=\bar{y}+r\frac{S_{Y}}{S_{X}}(x-\bar{x})$$

Propriétés : 

- La droite des moindres carrés passe par $(\bar{x} \ ; \ \bar{y})$, qui est le centre de gravité du nuage.       
- Le signe de sa pente est le même que celui de $S_{XY}$.

```{r Representation graphique du modele}
ggplot(df, aes(x=Superficie, y=Prix))+ 
  geom_point()+
  geom_smooth(method=lm, se=T)
```

### Contribution de chaque observation :

**Rappel :** la pente ($b$) de la droite des moindres carrés est calculée sur les observations $(x_{i} \ ; \ y_{i})$.

$$b= \frac{\sum_{i=1}^n (x_i-\bar x)(y_i-\bar y)}
       {\sum_{i=1}^n (x_i-\bar x)^2}
  =\sum^{n}_{i=1}\left[\frac{(x_{i}-\bar{x})^{2}}
  {\Sigma^{2}_{i=1}(x_{i}-\bar{x})^{2}}\right]\frac{(y_{i}-\bar{y})}{(x_{i}-\bar{x})}$$

Cette pente peut être interprétée comme **la moyenne pondérée des pentes des droites passant par le barycentre $(\bar{x}, \bar{y})$ du nuage et chaque point $(x_{i}, y_{i})$ individuellement.**

#### Levier : 

$$h_{i}=\frac{1}{n}+\frac{(x_{i}-\bar{x})^{2}}
  {\Sigma^{2}_{i=1}(x_{i}-\bar{x})^{2}}$$

Propriétés : 

- $h_{i}\in[\frac{1}{n}\ ;\ 1]$       
- $h_{i}=1 \ \Rightarrow \ \{\forall \ i \neq j \ ; \ x_{i}=x_{j}\}$      
- $h_{i}=\frac{1}{n} \ \Rightarrow \ \{\forall \ (j \neq k) \neq i\ ;\ x_{j}=x_{k}\}$
- Généralement, un levier suppérieur à $4/n$ est considéré comme important

**Attention :** le levier décrit un **potentiel**. Il ne décrit pas **directement** l'influence de l'observation. 

```{r Recherche de valeurs atypiques}
autoplot(
  modele,
  which = 3,
  ncol = 1,
  label.size = 3,
  label.hjust = -0.8,
  label.n = 6
) +
  aes(.fitted, .stdresid) + labs(x = "fitted values", y = "standardized residuals") +
  ggtitle("")

# Graphique basique des leviers :
plot(hatvalues(modele))

# Version amélioree du graphique des leviers :
ggplot(modele, aes(seq_along(.hat), .hat)) + geom_col(width = 0.1, colour = "blue") +
  labs(x = "Observation", y = "Leverage") + geom_text(
    label = rownames(df),
    check_overlap = T,
    vjust = -0.8,
    size = 3
  ) + geom_hline(yintercept = 4 / nrow(df),
                 colour = "red")

# Graphique pour la distznce de cook :  
ols_plot_cooksd_chart(modele)
 

ols_plot_dfbetas(modele)
ols_plot_dffits(modele)
ols_plot_resid_lev(modele)
```

### Résidus et qualité d'ajustement : 

```{r Etude des residus}
autoplot(
  modele,
  which = 1,
  ncol = 1,
  label.size = 3,
  label.hjust = -0.8,
  label.n = 6
)

autoplot(
  modele,
  which = 2,
  ncol = 1,
  label.size = 3,
  label.hjust = -0.8,
  label.n = 6
)
```

#### Décomposition de la variation totale :

$$SCT=\sum^{n}_{i=1}(y_{i}-\bar{y}_{i})^{2}=\sum^{n}_{i=1}(y_{i}-\hat{y}_{i}+\hat{y}_{i}-\bar{y}_{i})^{2}= \\ 
\sum^{n}_{i=1}(y_{i}-\hat{y}_{i})^{2}-2\sum^{n}_{i=1}(y_{i}-\hat{y}_{i})(\hat{y}_{i}-\bar{y}_{i})+\sum^{n}_{i=1}(y_{i}(\hat{y}_{i}-\bar{y}_{i})^{2} \\
$$

Avec :

$$
2\sum^{n}_{i=1}(y_{i}-\hat{y}_{i})(\hat{y}_{i}-\bar{y}_{i})=2\sum^{n}_{i=1}(y_{i}-\hat{y}_{i})\hat{y}_{i}-2\sum^{n}_{i=1}(y_{i}-\hat{y}_{i})\bar{y}_{i}
$$

Où :

$$
=2\sum^{n}_{i=1}(y_{i}-\hat{y}_{i})-2\bar{y}_{i}\sum^{n}_{i=1}(y_{i}-\hat{y}_{i})=2\sum^{n}_{i=1}(y_{i}-\hat{y}_{i})-2\bar{y}_{i}\sum^{n}_{i=1}e_{i}=2\sum^{n}_{i=1}(y_{i}-\hat{y}_{i}) \\
=2\sum^{n}_{i=1}e_{i}(a+bx_{i})=2a\sum^{n}_{i=1}e_{i}+b2\sum^{n}_{i=1}e_{i}x_{i}=0
$$

$$
SCT=\sum^{n}_{i=1}(y_{i}-\bar{y}_{i})^{2}=\sum^{n}_{i=1}(y_{i}-\hat{y}_{i})^{2}+\sum^{n}_{i=1}(\hat{y}_{i}-\bar{y}_{i})^{2}
$$

$$
SCT=SCE+SCR
$$

Avec : 

$$
SCE=\sum^{n}_{i=1}(\hat{y}_{i}-\bar{y}_{i})^{2}
$$

$$
SCR=\sum^{n}_{i=1}(y_{i}-\hat{y}_{i})^{2}$$
$$

```{r Affichage des parametres du modele}
anova(modele)
predict(modele,interval="confidence",level=0.95)
```

#### Coefficient de détermination : 

$$
R^{2}=\frac{SCE}{SCT}=\frac{\sum^{n}_{i=1}(\hat{y}_{i}-\bar{y}_{i})^{2}}{\sum^{n}_{i=1}(y_{i}-\hat{y}_{i})^{2}}
$$

Propriétés : 

- $R^{2}\in [\ 0 \ ; \ 1\ ]$
- $R^{2}=1 \ \Rightarrow$ liaison parfaite (les observations forment une droite).      
- $R^{2}=0 \ \Rightarrow$ aucune liaison linéaire.
- On démontre que $R^{2}$ est égal à la pente des moindres carrés.

### Propriétés des écarts résiduels :

Soit l'**écart résiduel** : $e_i = y_i - \hat{y}_i= y_i - \bigl[\bar{y} + b (x_i - \bar{x})\bigr]$

$$
\sum_{i=1}^{n}e_{i}=\sum_{i=1}^{n}(y_{i}-\hat{y}_{i})=\sum_{i=1}^{n}y_{i}-n\bar{y}-b\sum_{i=1}^{n}(x_{i}-\bar{x})=n\bar{y} - n\bar{y}=0
$$

Propriétés : 

- $E[e_{i}]=0$       
- La **variance empirique** des $e_{i}$ avec : 

$$S^{2}_{Y|X}=\frac{1}{n}\sum^{n}_{i=1}e_{i}^{2}=\frac{1}{n}\left[ \sum^{n}_{i=1}(y_{i}-\bar{y}_{i})^{2}-2b\sum_{i=1}^{n}(x_{i}-\bar{x})(y_{i}-\bar{y})+b^{2}\sum_{i=1}^{n}(x_{i}-\bar{x})^{2}\right] \\
=S^{2}_{Y}-2bS_{XY}+b^{2}S^{2}_{X}
$$

**Rappel :** 

- $b=r\frac{S_{Y}}{S_{X}}$       
- $r=\frac{S_{XY}}{S_{X}S_{Y}}$      

$$\Rightarrow \ S^{2}_{Y|X}=S^{2}_{Y}-2r\frac{S_{Y}}{S_{X}}rS_{X}S_{Y}+r^{2}S^{2}_{Y}=(1-r^{2})S^{2}_{Y}$$

### Coefficient de corrélation linéaire : 

$$\rho=\frac{Cov(X,Y)}{\sigma_{X}\sigma_{Y}}$$

Propriétés : 

- $\rho \in [-1 \ ; \ 1]$       
- Donne une indication sur le sens et l'intensité de la relation.       
- Généralement, $\rho$ mesure la qualité de l'approximation par une fonction linéaire.

### Dépendance entre deux variables : 

On sait que :

- $r=0 \ \nRightarrow\ Y \perp X$       
- $Y \perp X \ \Rightarrow\ r=0$       
- $Y \perp X \ \Leftrightarrow\ P(Y=y)=P(Y=y \ | \ X=x)$       

Pour être moins restrictif, on travail avec $E[\ Y|X\ ]$ afin de savoir si $Y$ dépend **en moyenne** de $X$.

#### Rapport de corrélation : 

- si $Y$ dépend en moyenne de $X$, alors il existe un **lien fonctionnel** du type : $Y=f(X)$.      
- On démontre que : $E[E[\ Y|X\ ]]=Y$       
- On démontre que : $V[\ Y|X\ ]=V[E[\ Y|X\ ]]+E[V[\ Y|X\ ]]$       

On déduit de ces observations le **rapport de corrélation** :

$$\eta^{2}_{Y|X}=\frac{V[E[\ Y|X \ ]]}{V[\ Y|X \ ]}$$

Propriétés : 

- $\eta \in [0 \ ; \ 1]$       
- $\eta = 0 \quad \Rightarrow$ indépendance en moyenne ; $E[\ Y|X \ ]$ est *certainement* constante.

- $\eta = 1 \ \Rightarrow E[V[\ Y|X \ ]]$ et $E[\ Y|X \ ]$ sont *certainement* nulles.         
- $E[\ Y|X \ ]=\alpha+X\beta$ (ex. $(X \space ; \space Y)$ Gaussiens) $\Rightarrow \space \eta_{Y | X}^{2}=\rho_{Y | X}^{2}$

### Modèle théorique :

Soit $X$ et $Y$, formant des couples d'observations tels que : $(X \ ; \ Y)=(X_{1} \ ; \ Y_{1}),...,(X_{n} \ ; \ Y_{n})$

On cherche $f$ telle que : $f(X) \approx Y \quad \Leftrightarrow \ E[\ (Y-f(X))^{2}\ ] \quad  \rightarrow 0$

On démontre que $E[\ (Y − f(X))^{2}]$ est minimale pour $f(X)=E[ \ Y |X \ ]$.

Soit la **fonction de régression** $f$ : 

$$f : X \rightarrow E[ \ Y |X=x \ ]$$

Alors : 

$$Y = E[ \ Y |X \ ] + \epsilon$$

Avec $\epsilon$, le résidu aléatoire tel que : 

- $E[\ \epsilon \ ]=0$ (conséquence de $E(Y ) = E [E(Y |X)]$)      
- $\epsilon \perp \{X \ ; \ E [Y |X] \}$       
- $V[\ \epsilon\ ]=V[\ Y-E[\ Y|X\ ] \ ]=(1-\eta^{2}_{Y|X})V[\ Y\ ]$
- $\epsilon \sim \mathcal{N}(0 \ ; \ \sigma^{2}) \quad \Leftrightarrow \quad Y_{i} \sim \mathcal{N}(\alpha+\beta x_{i} \ ; \ \sigma^{2})$ 


#### Régression linéaire : 

$$E[\ Y|X \ ]=\alpha+X\beta \ \ \Longrightarrow \ \ Y= \alpha+X\beta + \epsilon$$

Egalement : 

$$E[\ Y \ ] = E [E(Y |X)] \ \ \Longrightarrow \ \ E[\ Y \ ] = E[\ α + βX \ ] = α + βE[\ X \ ]$$

Donc la droite passe par le point $(\bar{X} \ ; \ \bar{Y})$. 

On a également : 

$$Y − E[ \ Y \ ] = α + βX + ϵ − [ \ α + βE(X \ ]] = β [ \ X − E(X) \ ] + \epsilon$$
$$E[\ (Y-E[\ Y\ ])(X-E[\ X\ ])\ ]=\beta E[\ (X-E[\ X\ ])^{2}\ ]+E[\ \epsilon(X-E[\ X\ ])\ ]$$

Sachant que : 

- $V[\ X \ ]=E[\ (X - E [ \ X\ ])^{2}\ ]$      
- $Cov(X \ ; \ Y)=E[\ (Y-E[\ Y\ ])(X-E[\ X\ ])\ ]$      
- $Cov(\epsilon \ ; \ X)=E[\ (\epsilon -E[\ \epsilon\ ])(X-E[\ X\ ])\ ]$      
- $E[\ \epsilon \ ]=0$

$$Cov(X \ ; \ Y)=\beta \ V[\ X\ ] + Cov(\epsilon \ ; \ X)$$

Alors même que $\epsilon$ n'est pas corrélée à $X$, soit $Cov(\epsilon \ ; \ X)=0$ et donc : 

$$Cov(X \ ; \ Y)=\beta \ V[\ X\ ]$$
$$\beta= \frac{Cov(X \ ; \ Y)}{ V[\ X\ ]}$$

Egalement : $\rho=\frac{Cov(X \ ; \ Y)}{\sigma_{X}\sigma_{Y}} \ \Rightarrow \ Cov(X \ ; \ Y) = \rho \sigma_{X}\sigma_{Y}$

$$\Rightarrow \ \ \beta=\frac{Cov(X \ ; \ Y)}{\sigma^{2}_{X}}=\frac{\rho \sigma_{X}\sigma_{Y}}{\sigma^{2}_{X}} \ \Rightarrow \ \beta=\rho\frac{\sigma_{Y}}{\sigma_{X}}$$

**Rappel :** L'équation de la droite de régression est donnée par :

$$Y-E[\ Y \ ]= \rho \frac{\sigma_{Y}}{\sigma_{X}}\left(X-E[\ X \ ]\right)+\epsilon$$
$$\Rightarrow \ Y= \left(E[\ Y \ ] - \rho \frac{\sigma_{Y}}{\sigma_{X}}E[\ X \ ]\right)+ \rho \frac{\sigma_{Y}}{\sigma_{X}}X + \epsilon$$

Si la régression est linéaire, on retrouve le résultat précédent ($\eta^{2}_{Y|X}=\rho^{2}$) :

$$
V[\ Y\ ]=V\left[\ E[\ Y \ ]+ \rho \frac{\sigma_{Y}}{\sigma_{X}}\left(X-E[\ X \ ]\right)+\epsilon\ \right]
$$

Avec $X$ et $\epsilon$ non corrélées :

$$
V[\ Y\ ]=V\left[\ \rho \frac{\sigma_{Y}}{\sigma_{X}}\left(X-E[\ X \ ]\right)+\epsilon\ \right] \\
=V\left[\ \rho \frac{\sigma_{Y}}{\sigma_{X}}\left(X-E[\ X \ ]\right) \right]+V\left[\ \epsilon\ \right]=\rho^{2}\frac{\sigma^{2}_{Y}}{\sigma^{2}_{X}}V\left[\ X \right]-V \left[\ E[\ X \ ]\right] +V\left[\ \epsilon\ \right]=\rho^{2}V\left[\ Y \right]+V\left[\ \epsilon\ \right]\\
V[\ \epsilon\ ]=(1-\rho^{2}) \ V[\ Y\ ]
$$

Or on a vu que : $V[\ \epsilon\ ]=(1-\eta^{2}_{Y|X}) \ V[\ Y\ ]$

Donc : 

$$\rho^{2}=\eta^{2}_{Y|X}$$  

### Modèle linéaire :

```{r Creation du modele}
modele=lm(Prix~Superficie, data = df)
```

*Les résultats obtenus jusqu'ici s'appliquent même si X n'est pas aléatoire mais contrôlée par l'expérimentateur*.

$$
\forall i : Y_{i}=\alpha + \beta x_{i} +\epsilon_{i}
$$

On parle de **modèle linéaire** et non de régression linéaire. 

Propriétés de la variable aléatoire $\epsilon$ :

- $\epsilon \sim \mathcal{N}(0 \ ; \ \sigma^{2})$
- $\forall i : V[\ \epsilon \ ]=\sigma^{2}$ (homoscédasticité)         
- $\epsilon_{i} \perp \epsilon_{j} \ , \forall i \neq j \ \Rightarrow \ Cov(\epsilon_{i},\epsilon_{j})=0$

Implications de l'hypothèse de normalité des résidus :

- $Y|X=x \sim \mathcal{N}(\alpha+\beta x \ ; \ \sigma^{2})$      
- $\forall x_{i} \ ; \hat Y_{i} \sim \mathcal{N}\left(\alpha+\beta x_{i} \ ; \ \sigma^{2}\left(\frac{1}{n}+\frac{(x_{i}-\bar{x})^{2}}{nS^{2}_{X}}\right)\right)$          
- $\forall x_{i} \ ; A \sim \mathcal{N}\left(\alpha \ ; \ \sigma^{2}\left(\frac{1}{n}+\frac{\bar{x}^{2}}{nS^{2}_{X}}\right)\right)$                
- $\forall x_{i} \ ; B \sim \mathcal{N}\left(\beta \ ; \ \frac{\sigma^{2}}{nS^{2}_{X}}\right)$  

```{r Evaluation des hypotheses du modele}
# Recherhe d'autocorrélation :
DurbinWatsonTest(modele, alternative = 'two.sided')

# Evaluation de l'homogeneite : 
bptest(modele)

# Test de la normalite :
shapiro.test(modele$residuals)
```
### Aspects inférentiels de la régression linéaire :

En partant du modèle :

$$Y_{i}=\alpha + \beta x_{i} +\epsilon_{i}$$

On effectue un développement par la méthode des moindres carrés, comparable à celui déja effectué, afin de calculer les paramètres en tenant compte du caractère aléatoire de $Y$ :

$$B=\frac{\sum^{n}_{i=1} (Y_{i}-\bar{Y})(x_{i}-\bar{x})}{\sum^{n}_{i=1}(x_{i}-\bar{x})^{2}}$$
$$A=\bar{Y}-B\bar{x}$$

#### Propriétés de $A$ et $B$ :

```{r Affichage des coefficients}
summary(modele)
```

- **Linéarité** car ils sont combinaisons linéaires des $Y_{i}$ :

$$\delta_{i}=\sum^{n}_{i=1}\frac{(x_{i}-\bar{x})}{\sum^{n}_{i=1}(x_{i}-\bar{x})^{2}}Y_{i.} \ \qquad \Longrightarrow \qquad B=\Sigma^{n}_{i=1}\delta_{i}Y_{i} \qquad et \qquad A=\Sigma^{n}_{i=1}\left(\frac{1}{n}\bar{x}\delta_{i}\right)Y_{i}$$

- Absence de **biais** : 

$$E[\ B \]=E[\frac{\sum_{i=1}^{n} (Y_i-\bar Y)(x_i-\bar x)}{\sum_{i=1}^{n}(x_i-\bar x)^{2}}\]=  \frac{\sum^{n}_{i=1} E[Y_{i}-\bar Y](x_{i}-\bar x )}{\sum_{i=1}^{n}(x_{i}-\bar{x})^{2}}$$

On sait que :

$$E[\ Y_{i} \]=E[\ \alpha+\beta x_{i}+\epsilon_{i} \]=\alpha+\beta x_{i}$$

$$E[\ \bar{Y} \ ]=\frac{1}{n}\Sigma_{i}E[\ Y_{i} \ ]=\frac{1}{n}\Sigma_{i}(\alpha+\beta x_{i})=\alpha+\beta \bar{x}$$

$$\Rightarrow E[\ B \ ]=  E\left[\ \frac{\sum^{n}_{i=1} \beta (x_{i}-\bar{x})(x_{i}-\bar{x})}{\sum^{n}_{i=1}(x_{i}-\bar{x})^{2}} \ \right]=\beta$$

$$
E[\ A\ ]=E[\ \bar{Y}-\beta\bar{x}\ ]=E[\ Y\ ]-\beta x=(\alpha+\beta x+E[\ \epsilon\ ])-\beta x=\alpha
$$

On démontre que : 

$$
V[\ B\ ]=\frac{\sigma^{2}}{nS^{2}_{X}}
$$

$$
V[\ A\ ]=\frac{\sigma^{2}}{n}\left(1+\frac{\bar{x}}{S^{2}_{X}}\right)
$$

Où $\sigma^{2}$ est la variance des résidus : $S_{X}^{2}=\frac{1}{n}\Sigma_{i}(x_{i}-\bar{x})^{2}$

- **non indépendance** : $A \not\perp B$

On démontre que : 

$$Cov(A \ ; \ B)=\frac{-\bar{x}\sigma^{2}}{\Sigma^{n}_{i=1}(x_{i}-\bar{x})^{2}}$$

$$\rho(A \ ; \ B)=\frac{-\bar{x}}{\sqrt{\Sigma^{n}_{i=1}(x_{i}^{2})/n}}$$

- $A$ et $B$ sont des estimateurs de **variance minimale** (théorème de *Gauss-Markov*).

### Estimation par maximum de vraissemblance : 

**Rappel :**

- $\epsilon \sim \mathcal{N}(0 \ ; \ \sigma^{2})$

Soit un échantillon $(Y_{1},...,Y_{i},...,Y_{n})$ avec $Y_{i} \sim \mathcal{N}(\alpha + \beta x_{i} \ ; \ \sigma^{2})$.

Soit une fonction de densité de $Y_{i}$ :

$$f(Y_{i})=\frac{1}{\sqrt{2\pi\sigma^{2}}}e^{-\frac{(Y_{i}-\alpha-\beta x_{i})^{2}}{2\sigma^{2}}}$$

Avec la fonction de vraissemblance : 

$$L(\alpha, \beta,\sigma^{2}, Y_{1},...,Y_{n})=
\prod_{i=1}^{n} \frac{1}{\sqrt{2\pi\sigma^{2}}}e^{-\frac{(Y_{i}-\alpha-\beta x_{i})^{2}}{2\sigma^{2}}}$$

$$L(\alpha, \beta,\sigma^{2}, Y_{1},...,Y_{n})=
\left(\frac{1}{\sqrt{2\pi\sigma^{2}}}\right)^{n}e^{-\frac{1}{2\sigma^{2}}\sum^{n}_{i=1}(Y_{i}-\alpha-\beta x_{i})^{2}}$$

$$\mathcal{l}(\alpha, \beta,\sigma^{2}, Y_{1},...,Y_{n})=-ln((2\pi\sigma^{2})^{\frac{n}{2}})-\frac{1}{2\sigma^{2}}\sum^{n}_{i=1}(Y_{i}-\alpha-\beta x_{i})^{2}$$

$$\mathcal{l}(\alpha, \beta,\sigma^{2}, Y_{1},...,Y_{n})=-\frac{n}{2}ln(2\pi)-\frac{n}{2}ln(\sigma^{2})-\frac{1}{2\sigma^{2}}\sum^{n}_{i=1}(Y_{i}-\alpha-\beta x_{i})^{2} \\
$$

On recherche un extrêmum :

$$
\frac{\partial ln(L)}{\partial \alpha}=\frac{1}{\sigma^{2}}\sum^{n}_{i=1}(Y_{i}-\alpha-\beta x_{i})$$

$$\frac{\partial ln(L)}{\partial \beta}=\frac{1}{\sigma^{2}}\sum^{n}_{i=1}(Y_{i}x_{i}-\alpha x_{i}-\beta x_{i}^{2})$$

$$\frac{\partial ln(L)}{\partial \sigma^{2}}=-\frac{n}{2\sigma^{2}}-\frac{1}{2\sigma^{4}}\sum^{n}_{i=1}(Y_{i}x_{i}-\alpha x_{i}-\beta x_{i}^{2})
$$

Donc :

$$\sum^{n}_{i=1}(Y_{i}-\alpha-\beta x_{i}) = 0 \quad \Rightarrow \quad A=\bar{Y}-B\bar{X} $$
$$\sum^{n}_{i=1}(Y_{i}x_{i}-\alpha x_{i}-\beta x_{i}^{2}) = 0 \quad \Rightarrow \quad B= \frac{\sum^{n}_{i=1} (Y_{i}-\bar{Y})(x_{i}-\bar{x})}{\sum^{n}_{i=1}(x_{i}-\bar{x})^{2}} $$
$$-\frac{n}{2\sigma^{2}}-\frac{1}{2\sigma^{4}}\sum^{n}_{i=1}(Y_{i}x_{i}-\alpha x_{i}-\beta x_{i}^{2}) = 0 \quad \Rightarrow \quad S^{2}_{n}=\frac{1}{n}\sum^{n}_{i=1}(Y_{i}-\alpha-\beta x_{i})^{2} $$

Remarque : 

- Les estimateurs $\alpha$ et $\beta$ coïncident avec les estimateurs des moindres carrés.     
- L'EMV de $\sigma^{2}$ est **biaisé**. l'estimaeur sans biais de $\sigma^{2}$ est :

$$S^{2}_{n-2}=\frac{1}{n-2}\sum^{n}_{i=1}(Y_{i}-\hat{Y_{i}})^{2}=\frac{1}{n-2}\sum^{n}_{i=1}(Y_{i}-A-Bx_{i})^{2}$$

Son estimation est donc : 

$$
\hat{\sigma}^{2}=\frac{1}{n-2}\sum^{n}_{i=1}(y_{i}-\hat{y_{i}})^{2}=\frac{1}{n-2}\sum^{n}_{i=1}e_{i}^{2}
$$

### Lois standardisées de $A$ et $B$ :

$$\frac{A-E[\ A\ ]}{V[\ A \ ]}=\frac{A-\alpha}{\sqrt{\frac{\sigma^{2}}{n}\left(1+\frac{\bar{x}^{2}}{S^{2}_{X}} \right)}}=\frac{A-\alpha}{\sqrt{\sigma^{2}\left(\frac{1}{n}+\frac{\bar{x}^{2}}{\sum^{n}_{i=1}(x_{i}-\hat{x})^{2}} \right)}} \sim \mathcal{N}(0 \ ; \ 1)$$
$$\frac{B-E[\ B\ ]}{V[\ B \ ]}=\frac{B-\beta}{\sqrt{\frac{\sigma^{2}}{nS^{2}_{X}}}}=\frac{B-\beta}{\sqrt{\frac{\sigma^{2}}{\sum^{n}_{i=1}(x_{i}-\hat{x})^{2}}}} \sim \mathcal{N}(0 \ ; \ 1)
$$

$$
\frac{(n-2)S^{2}_{n-2}}{\sigma^{2}}=\frac{\sum^{n}_{i=1}(Y_{i}-A-Bx_{i})^{2}}{\sigma^{2}}=\frac{SCR}{\sigma^{2}} \sim \chi^{2}_{n-2}
$$

On démontre que $A \perp B$

A partir de ces observations, on trouve les lois *studentisées* de $A$ et $B$ :

$$
\frac{A-\alpha}{S_{A}} \sim \mathcal{T}_{n-2}
$$

Car : $S_{A}=S_{n-2}\sqrt{\frac{1}{n}+\frac{\bar{x}^{2}}{\sum_{i=1}^{n}(x_{i}-\hat{x})^{2}}}$

$$
\frac{A-\alpha}{S_{A}}=\frac{A-\alpha}{\sqrt{\frac{(n-2)S^{2}_{n-2}}{\sigma^{2}(n-2)}}\sigma\sqrt{\frac{1}{n}+\frac{\bar{x}^{2}}{\sum^{n}_{i=1}(x_{i}-\hat{x})^{2}}}}
$$

Est de la forme : 

$$
\frac{A-\alpha}{\frac{\sqrt{\frac{\sigma^{2}}{n}\left(1+\frac{\bar{x}^{2}}{\sum^{n}_{i=1}(x_{i}-\hat{x})^{2}}\right)}}{\sqrt{\frac{\frac{(n-2)S^{2}_{n-2}}{\sigma^{2}}}{(n-2)}}}}
$$

Où le numérateur suit une loi normale ( $\mathcal{N} (0 \space ; \space 1)$ ), et le dénominateur est la racine d'une loi du $\chi^{2}$ divisée par son degrés de liberté $(n-2)$.

De la même façon, si on note : $S_{B}=\frac{S_{n-2}}{\sqrt{\sum_{i=1}^{n}(x_{i}-\bar{x})^{2}}}$

Alors : 

$$
\frac{B-\beta}{S_{B}} \sim \mathcal{T}_{n-2}
$$

**Remarque :** ces variables dépendent du paramètre alors que leurs lois en sont indépendantes (quantités **pivotales**). On peut donc les utiliser pour construire des tests statistiques. 

### Intervalles de confiance et de prédiction :

#### Détermination d'Intervalles de confiance pour $\alpha$ et $\beta$ :

**Rappel :** 

- $\frac{A-\alpha}{S_{A}} \sim \mathcal{T}_{n-2}$          
- $\frac{B-\beta}{S_{B}} \sim \mathcal{T}_{n-2}$

Soit $t_{n-2 \ ; 1-\gamma/2}$ , le fractile d'ordre $(1-\gamma)$ d'une loi de Student à $(n-2)$ degrés de liberté : 

$$P\left(\ -t_{n-2 \ ; 1-\gamma/2}\leq\frac{A-\alpha}{S_{A}}\leq t_{n-2 \ ; 1-\gamma/2}\right)=1-\gamma$$
$$P\left(\ -t_{n-2 \ ; 1-\gamma/2}\times S_{A}\leq A-\alpha\leq t_{n-2 \ ; 1-\gamma/2}\times S_{A}\right)=1-\gamma$$
$$P\left(\ -A -t_{n-2 \ ; 1-\gamma/2}\times S_{A}\leq -\alpha\leq -A +t_{n-2 \ ; 1-\gamma/2}\times S_{A}\right)=1-\gamma$$
$$P\left(\ A -t_{n-2 \ ; 1-\gamma/2}\times S_{A}\leq \alpha\leq A +t_{n-2 \ ; 1-\gamma/2}\times S_{A}\right)=1-\gamma
$$

Donc : 

$$
IC_{n-2 \ ; \ (1-\gamma/2)}=\left[A -t_{n-2 \ ; 1-\gamma/2}\times \hat{\sigma}_{A} \quad ; \quad A +t_{n-2 \ ; 1-\gamma/2}\times \hat{\sigma}_{A} \right]
$$

et : 

$$P\left(\ -t_{n-2 \ ; 1-\gamma/2}\leq\frac{B-\beta}{S_{B}}\leq t_{n-2 \ ; 1-\gamma/2}\right)=1-\gamma$$
$$P\left(\ B -t_{n-2 \ ; 1-\gamma/2}\times S_{B}\leq \beta\leq B +t_{n-2 \ ; 1-\gamma/2}\times S_{B}\right)=1-\gamma$$

$$IC_{n-2 \ ; \ (1-\gamma/2)}=\left[B -t_{n-2 \ ; 1-\gamma/2}\times \hat{\sigma}_{B} \quad ; \quad B +t_{n-2 \ ; 1-\gamma/2}\times \hat{\sigma}_{B} \right]
$$

#### Détermination d'Intervalles de confiance pour la droite de régression :

Soit : $y_{i} \ : \ E[\ Y | X= x_{i }\ ]=\alpha+\beta x_{i}$         
dont l'estimateur est : $A+Bx_{i}$

On a :

$$E[\ A+Bx_{i}\ ]=E[\ A\ ]+E[\ B\ ]x_{i}=\alpha+\beta x_{i}$$

$$V[\ A+Bx_{i}\ ]=V[\ A\ ]+V[\ B\ ]x_{i}^{2}+2x_{i}\ Cov(A\ ;B)$$
$$=\sigma^{2}\left[\frac{1}{n}+\frac{\bar{x}^{2}}{\sum^{n}_{j=1}(x_{j}-\bar{x})^{2}}\right]+x_{i}^{2}\frac{\sigma^{2}}{\sum^{n}_{j=1}(x_{j}\bar{x})^{2}}-2x_{i}\frac{\bar{x}\sigma^{2}}{\sum^{n}_{j=1}(x_{j}-\bar{x})^{2}}$$
$$=\sigma^{2}\left[\frac{1}{n}+\frac{(x_{i}-\bar{x})^{2}}{\sum^{n}_{j=1}(x_{j}-\bar{x})^{2}}\right]=\sigma^{2}h_{i}$$

La quantité $A+B x_{i}$ pouvant être exprimée comme une **combinaison linéaire de variables aléatoires gaussiènnes**, elle suit une loi normale : 

$$
A+B x_{i} \sim \mathcal{N}(\alpha+\beta x_{i} \ ; \ \sigma \sqrt{h_{i}})
$$

Donc :

$$
\frac{A+Bx_{i}-(\alpha+\beta x_{i})}{\sigma\sqrt{\frac{1}{n}+\frac{(x_{i}-\bar{x})^{2}}{\sum^{n}_{j=1}(x_{j}-\bar{x})^{2}}}} \sim \mathcal{N}(0 \ ; \ 1)
$$

Soit en utilisant l'estimateur de $\sigma$ :

$$
\frac{A+Bx_{i}-(\alpha+\beta x_{i})}{S_{n-2}\sqrt{\frac{1}{n}+\frac{(x_{i}-\bar{x})^{2}}{\sum^{n}_{j=1}(x_{j}-\bar{x})^{2}}}} \sim \mathcal{T}_{n-2}
$$

Avec :

$$
IC_{1-\gamma} \quad E[\ Y|x_{i}\ ]=\left[\ \hat{y}_{i} \ -\mathcal{T}_{n-2 \ ; \ (1-\gamma/2)} \times \hat{\sigma}\sqrt{h_{i}} \quad ;\quad \hat{y}_{i} \ +\mathcal{T}_{n-2 \ ; \ (1-\gamma/2)} \times \hat{\sigma}\sqrt{h_{i}} \ \right]
$$

Remarque : l'éloignement de $x_{i}$ par rapport à $\bar{x}$ tend à élargir l'IC.

#### Intervalle de prévision : 

Soit la valeur non observée $x_{0}$ et le modèle de prévision suivant :

$$
y_{0}=\alpha+\beta x_{0} +\epsilon_{0}
$$

**Rappel :** $y_{0}$ et $\epsilon_{0}$ sont des variables aléatoires.

Soit $y^{*}$, la prévision naturelle de $Y$ :

$$
y^{*}_{0}=a+bx_{0}
$$

On cherche ($Y_0-\hat Y_{0}$), le résidus d'une nouvelle observation pour laquelle on dispose de $x_{0}$. On a vu qu'en $x_{0}$ :

$$
\hat{Y}=A+Bx_{0} \sim \mathcal{N}(\alpha+\beta x_{0} \ ;\ \sigma^{2}\sqrt{h_{0}})
$$

Alors que :

$$
Y|X=x_{0} \sim \mathcal{N}(\alpha+\beta x_{0} \ ; \ \sigma)
$$

On remarque qu'alors que $Y_{0}$ ne dépend que de $x_{0}$, la valeur de $\hat{Y}_{0}$ dépend quand à elle des **valeurs déja observées**. 

Donc : 

$$
V[\ (Y-\hat{Y})\ ]=V[\ Y \ ]-V[\ \hat{Y}\ ]+2 \ Cov(Y\ ; \ \hat{Y})=V[\ Y \ ]-V[\ \hat{Y} = \sigma^{2}+\sigma^{2}h_{0}=\sigma^{2}(1+h_{0}
$$

$$(Y_{0}-\hat{Y}_{0}) \sim \mathcal{N}(0 \ ; \sigma\sqrt{1+h_{0}})$$

On standardise et on utilise l'estimateur de $\sigma^{2}$ :

$$\frac{Y_{0}-\hat{Y}^{\*}_{0}}{S_{n-2}\sqrt{1+h_{0}}} \sim \mathcal{T}_{n-2}$$

$$IC_{1-\gamma} \ (Y_{0}-\hat{Y^{\*}}_{0})=IC_{1-\gamma} \ [\ Y|x_{0}\ ]=\left[\ \hat{Y^{\*}}_{0}-\mathcal{T}_{n-2 \ ; \ (1-\gamma/2)} \times S_{n-2}\sqrt{1+h_{0}} \ ; \ \hat{Y^{\*}}_{0}+\mathcal{T}_{n-2 \ ; \ (1-\gamma/2)} \times S_{n-2}\sqrt{1+h_{0}}\ \right]$$

Remarque : On trouve la même tendance qu'avec l'intervalle de confiance de la droite, avec un élargissement proportionnel à l'éloignement vis à vis de $\bar{x}$.

```{r Interval de prediction}
PI = as.data.frame(cbind(Prix = df$Prix, Superficie = df$Superficie,
predict(modele,interval="prediction")))

ggplot(PI, aes(x=Superficie, y=Prix))+
geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
geom_line(aes(y=upr), color = "red", linetype = "dashed")+
geom_point()+geom_smooth(method=lm, se=T) + geom_text(label = row.names(df),
vjust = - 1, check_overlap = TRUE, size = 3)
```

### Tests d'hypothèse sur les paramètres du modèle :

$$
\frac{B-\beta}{S_{B}} \sim \mathcal{T}_{n-2} \quad SSi \quad \beta = \beta_{0}
$$

$h_{0} \ : \ \beta = \beta_{0}$          
$h_{1} \ : \ \beta \neq \beta_{0}$

Remarque : l'hypothèse alternative est **bilatérale**.

Statistique : 

$$
T_{obs}=\frac{b}{\hat{\sigma}_{B}}
$$

Egalement :        

$$
\frac{A-\alpha}{S_A} \sim \mathcal{T}_{n-2} \quad SSi \quad \alpha = \alpha_0
$$

$h_{0} \ : \ \alpha = \alpha_{0}$          
$h_{1} \ : \ \alpha \neq \alpha_{0}$

Remarque : l'hypothèse alternative est **bilatérale**.

Statistique : 

$$
T_{obs}=\frac{a}{\hat{\sigma}_{A}}
$$

#### Test de significativité globale du modèle :

$H_{0} \ : \ y_{i}=\alpha +\epsilon_{i}$ (modèle nul) est le vrai modèle.                
$H_{1} \ : \ y_{i}=\alpha +\epsilon_{i}$ '(modèle nul) n'est pas le vrai modèle.                

**Rappel :** $\frac{\sum_{i=1}^{n}(Y_{i}^{*}-\bar Y_{i})^{2}}{\sigma^{2}}=\frac{SCE}{\sigma^{2}} \quad \sim_{H_{0}} \quad \chi_{1}^{2}$

On démontre également que $SCE \perp SCR$

On fabrique donc la statistique suivante : 

$$
\frac{SCE}{SCR/(n-2)} \quad \sim_{H_{0}} \quad \mathcal{F}_{(1 \ ; \ n-2)}
$$

Ce test s'accompange généralement d'un **tableau d'analyse de variance** : 

| Source   | DDL     | Somme des carrés                              | Carrés moyens        | Statistique du test                     | P-value              |
|----------|---------|-----------------------------------------------|----------------------|------------------------------------------|----------------------|
| Modèle   | 1       | SCE = $\sum_{i=1}^n (\hat Y_i - \bar Y)^2$     | $\dfrac{\text{SCE}}{1}$ | $F = \dfrac{\text{SCE}}{\text{SCR}/(n-2)}$ | $P(F > F_{\text{obs}})$ |
| Résidu   | $n-2$   | SCR = $\sum_{i=1}^n (Y_i - \hat Y_i)^2$        | $\dfrac{\text{SCR}}{n-2}$ |                                          |                      |
| Totale   | $n-1$   | SCT = $\sum_{i=1}^n (Y_i - \bar Y)^2$          | $\dfrac{\text{SCT}}{n-1}$ |                                          |                      |

### Etude des résidus :

#### Différence entres erreurs aléatoires, résidus aléatoires et résidus observés :

- Erreur aléatoire : $\epsilon_{i}=Y_{i}-\alpha-\beta x_{i}$ (variable aléatoire)        
- Résidus aléatoire : $E_{i}=Y_{i}-A-Bx_{i}=Y_{i}-\hat{Y}_{i}$ (estimateur)
- Résidu observé : $e_i = y_i - \hat{y}_i$ (estimation de $\epsilon_i$)

Propriétés des résidus aléatoires : 

$$V[\ E_{i}\ ]=V[\ Y_{i}-[A+Bx_{i}]\ ]=V[\ Y_{i}\ ]+V[\ (A+Bx_{i})\ ]- 2 \ Cov(Y_{i} \ ; A+Bx_{i}\ )$$
$$= V[\ Y_{i}\ ]+V[\ (A+Bx_{i})\ ]- 2 \ V[\ \hat{Y_{i}} \ ] = \sigma^{2}+\sigma^{2}h_{i}-2\sigma^{2}h_{i}$$

$$\sigma^{2}_{E_{i}}= (1-h_{i}) \ \sigma^{2}$$

#### Résidus normalisés : 

Utilisés pour l'observation à partir d'une variance unitaire (indépendante des échelles de variables) :

- Standardisation interne (studentisés interne) :         

$$
\frac{e_{i}}{\hat{\sigma}\sqrt{1-h_{i}}}
$$

- Studentisés (studentisés externe) : Calcul de la variance sans utiliser l'observation $i$. Si on suppose vrai l'**hypothèse de normalité**. alors :         

$$
\frac{e_{i}}{\hat{\sigma}_{(-i)}\sqrt{1-h_{i}}} \quad  \sim \quad \mathcal{T}_{n-2}
$$

#### Graphe des résidus :

Permet de vérifier visuellement certaines hypothèses : linéarité, indépendance, homoscédasticité.

**Evaluation de l'hypothèse de normalité :** 

- graphique : QQ-plot        
- tests statistiques : Kolmogorov-Smirnov, Kramer-Von Mises, Anderson-Darling, Shapiro-Wilk.

**Evaluation de l'hypothèse d'indépendance et recherche de corrélations :** 

Test de Durbin-Watson : régression du type $\epsilon_{i+1}=\rho\epsilon_{i}+u_{i+1}$

- $H_{0} \ : \ \rho=0$                   
- $H_{1} \ : \ \rho\neq0$

$$
DW=\frac{\sum^{n}_{i=1}(e_{i+1}-e_{i})^{2}}{\sum^{n}_{i=1}e_{i}^{2}}
$$

Propriétés : 

- $DW\in [0 \ ; \ 4]$, avec une auto-corrélatino négative pour ($DW<2$) ou positive pour ($DW>2$).

### Recherche d'observations influentes : 

Critères de présélection : 

- Valeur atypique sur Y donne un résidu élevé.               
- valeur atypique sur X donne un levier important.                      

#### Distance de Cook :

$$
D_{i}=\sum^{n}_{i=1}\frac{(\hat{y}_{i}-\hat{y}_{(-i)j})^{2}}{2\hat{\sigma}^{2}}
$$

On peut l'écrire différement : 

$$
D_{i}=\frac{e_{i}^{2}}{2\hat{\sigma}^{2}}\times \frac{h_{i}}{(1-h_{i})^{2}}
$$

### Autres mesures à documenter :### 

Il existe d’autres mesures d’influence plus spécifiques :

— $Covratio_{i}$ mesure l’influence d’une observation sur la variabilité des estimations. $Covratio_{i}$ est calculé comme le rapport entre la variance des estimations sans et avec l’observation. On le considère important si $|Covratio_{i}-1| > 6/n$.

— $DFBETAS_{i}$ mesure l’influence d’une observation sur l’estimation des coefficients.
$DFBETAS_{i}$ est calculés comme différence normalisée entre la valeur de l’estimation
de $\alpha$ (ou $\beta$) avec et sans l’observation. On le considère important si $|DFBETAS| >\frac{2}{\sqrt{n}}$

— $DFFITS_{i}$, qui mesure l’influence d’une observation sur l’estimation de $E(Y |X = xi)$.
$DFFITS_{i}$ est calculés comme différence normalisée entre la valeur de l’estimation de
$E(Y |X = xi)$ avec et sans l’observation. On le considère important si $|DFFITS| >2\sqrt{2/n}$.


### Références :

[1] Gilbert Saporta. Probabilités, analyse des données et statistique. Editions Technip, 2006.

[2] Renée Veysseyre. Aide-mémoire-Statistique et probabilités pour les ingénieurs. Dunod, 2014. 34

[3] Giorgio Russolillo. STA102 : La régression Linéaire Simple (présentation ppt)

[4] Michel Lejeune. Statistique - La théorie et ses applications. Springer-Verlag. 35


```{r }
confint(modele)

modele_df = as.data.frame(cbind(
  df,
  predict(modele, interval = "confidence", level = 0.95)
))
```
























































































