Suppose we have $n$ samples from the following linear model:
$$
\begin{align}
y=x_1 \beta_1+x_2 \beta_2 + e,
\end{align}
$$
where $e_i$  i.i.d. comes from $N(0,\sigma^2)$; $x_1$, $x_2$ are centralized with mean $0$. 

I am **only interested in $\beta_1$**. 

1. While, we know that the least squares estimate $\widehat{\beta}_1^{(1)} =(x_1^\prime x_1)^{-1} x_1^\prime y$ is poor, because it doesn't consider the effect of $x_2$ on $y$ and could lead to large estimation variance. 
2. Let $x = (x_1, x_2)\in\mathbb{R}^{n\times 2}$. 
We have the least-squares estimate of $(\beta_1, \beta_2)^\prime$,
$$
\begin{align}
(\widehat{\beta}_1^{(2)} , \widehat{\beta}_2^{(2)})^\prime
& = (x^\prime x)^{-1} x^\prime y.
\end{align}
$$

3. Denote by $P_M$ the projection matrix corresponding to a design matrix $M$. 
    Denote by $I_n$ the identity matrix of order $n$. 
    Let
$$
\begin{align}
  \widehat{\beta}_1^{(3)} 
  & = \{x_1^\prime [I_n − x_2 (x_2^\prime x_2 ^{−1}) x_2^\prime ] x_1 \}^{−1} x_1^\prime [I_n−x_2 (x_2^\prime x_2 )^{−1} x_2^\prime ] y \\
  & = [x_1^\prime (I_n − P_{x_2} ) x_1 ]^{−1} x_1^\prime (I_n−P_{x_2} ) y \\
  & = (x_1^\prime H^\prime H X_1)^{-1} x_1^\prime H^\prime y \\
  & = \mathop{\arg\min}_{\beta_1} \| H (y - x_1\beta_1) \|_2^2 
  \end{align}
$$
where $H=(I_n - P_{x_2})$​ is also a projection matrix. 

$H x_1$: A projection of $x_1$ to the orthogonal complement space of column space of $x_2$​​. 

The column space of $x_2$​ is compressed, but the relationship between  $x_1$​ and  $x_2$​ is preserved.


$$
\begin{align}
\| y - x_1 \beta_1 \|_2^2 
& = \| (I_n − P_{x_2} ) (y - x_1\beta_1) + P_{x_2} (y - x_1\beta_1)\|_2^2 \\
& = \| (I_n − P_{x_2} ) (y - x_1\beta_1) \|_2^2 + \|P_{x_2} (y - x_1\beta_1)\|_2^2
\end{align}
$$

Let consider the linear equation problem (noiseless case): 
$$
\begin{align}
x^\prime x \beta & = x^\prime \mu \\
\begin{pmatrix}
x_1^\prime x_1 & x_1^\prime x_2 \\
x_2^\prime x_1 & x_2^\prime x_2 \\
\end{pmatrix}
\begin{pmatrix}\beta_1 \\ \beta_2 \end{pmatrix} & = \begin{pmatrix}x_1^\prime \mu \\ x_2^\prime \mu \end{pmatrix}
\label{eq:tmp11}
\end{align}
$$
By the second row in $\eqref{eq:tmp11}$, we have 
$$
\begin{align}
x_2^\prime x_1 \beta_1 + x_2^\prime x_2 \beta_2  & =  x_2^\prime \mu \\
\beta_2 & = (x_2^\prime x_2)^{-1} (x_2^\prime \mu - x_2^\prime x_1 \beta_1)
\end{align}
$$
We plug $\beta_2$ into the first row, in $\eqref{eq:tmp11}$: ​
$$
\begin{align}
x_1^\prime x_1 \beta_1 + x_1^\prime x_2 \beta_2  & =  x_1^\prime \mu \\
x_1^\prime x_1 \beta_1 + x_1^\prime x_2 (x_2^\prime x_2)^{-1} (x_2^\prime \mu - x_2^\prime x_1 \beta_1)  & =  x_1^\prime \mu  \\
x_1^\prime( I_n - P_{x_2})x_1\beta_1 &= x_1^\prime(I_n - P_{x_2}) \mu \\
\end{align}
$$
It leads to
$$
\begin{align}
\beta_1 &= (x_1^\prime( I_n - P_{x_2})x_1)^{-1} x_1^\prime(I_n - P_{x_2}) \mu .
\end{align}
$$

Replacing  $u$ by $y$, we have,
$$
\begin{align}
\widehat{\beta}_1^{(3)} 
& = (x_1^\prime( I_n - P_{x_2})x_1)^{-1} x_1^\prime(I_n - P_{x_2}) (\mu + e) \\
& = \beta_2 + (x_1^\prime( I_n - P_{x_2})x_1)^{-1} x_1^\prime(I_n - P_{x_2}) e .
\end{align}
$$

However, the LSE of $\beta_1$ is 
$$
\begin{align}
\widehat{\beta}_1^{(2)} 
& = \begin{pmatrix}1 & 0\end{pmatrix}
    \begin{pmatrix} x_1^\prime x_1 & x_1^\prime x_2 \\ x_2^\prime x_1 & x_2^\prime x_2 \\ \end{pmatrix}^{-1}
    \begin{pmatrix} x_1^\prime y \\ x_2^\prime y\\ \end{pmatrix} \\
& = \beta_1 + \begin{pmatrix}1 & 0\end{pmatrix}
    \begin{pmatrix} x_1^\prime x_1 & x_1^\prime x_2 \\ x_2^\prime x_1 & x_2^\prime x_2 \\ \end{pmatrix}^{-1}
    \begin{pmatrix} x_1^\prime e \\ x_2^\prime e \\ \end{pmatrix} \\
& = \beta_1 + (x_1^\prime x_1 x_2^\prime x_2 - x_1^\prime x_2 x_1^\prime x_2)^{-1} (x_2^\prime x_2 x_1^\prime - x_2^\prime x_1 x_1^\prime) e \\
& = \beta_1 + (x_1^\prime (x_1 x_2^\prime - x_2 x_1^\prime) x_2)^{-1} x_2^\prime(x_2  - x_1)x_1^\prime e \\
\end{align}
$$
Therefore, $\widehat{\beta}_1^{(3)}\neq\widehat{\beta}_1^{(2)}$​ although they are both unbiased estimate of $\beta_1$​​​​. Because their estimation errors are different. 

My **questions** are: 

1. Is $\widehat{\beta}_1^{(3)}$ also the least-squares estimate of $\beta_1$? 
    If not, what is $\widehat{\beta}_1^{(3)}$ and what does it mean? 

2. What are the differences between $\widehat{\beta}_1^{(1)}$​, $\widehat{\beta}_1^{(2)}$​ and $\widehat{\beta}_1^{(3)}$​?

   $\widehat{\beta}_1^{(2)} \neq \widehat{\beta}_1^{(3)}$​ in theory: 

   

   

---
I am very curious about $\widehat{\beta}_1^{(3)}$. But I was not able to find relevant materials for self-study. 

Although, I think $\widehat{\beta}_1^{(3)}$ is a weighted least squares estimate with weight $[I_n−x_2 (x_2^\prime x_2 )^{−1} x_2^\prime ]$. 
However, I don't know how to derive this [weighted least-squares estimate](https://en.wikipedia.org/wiki/Weighted_least_squares), what does the weight matrix mean, and when to use $\widehat{\beta}_1^{(2)}$ or $\widehat{\beta}_1^{(3)}$. 

Thanks. 