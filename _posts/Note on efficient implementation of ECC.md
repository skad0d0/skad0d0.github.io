---
title: 'Notes on Efficient implementation of ECC'
date: 2024-02-15
permalink: /posts/2024/02/blog-post-1/
tags:
  - ECC
---
## Sandy2x: New Record on Curve 25519
We define Curve25519 as a function that maps two 32-byte input string to a 32-byte output string. The function can be viewed as a x-coordinate only scalar multiplication on the curve 
$$ E_M : y^2 = x^3 + 486662x^2 + x $$
over $ \mathbb{F}_{2^{255}-19} $. One input string is an integer scalar $s$ i.e. the public string. The other input string is a 32-byte encoding of $x_P$ i.e. the x-coordinate of a point P on the curve.

Here we give a high level understanding of the function Curve25519. For a set of two Curve25519 users, they each have a 32-byte secret key and 32-byte public key, and they have a 32-byte **shared secret key** to authenticate and encrypt message between the two users. E.g. Alice holds a secrete key $a$, and Bob holds a secret key $b$, and the public string is $9$. Then Alice's public key is $Curve25519(a,9)$, and Bob's public key is $Curve25519(b,9)$. For the shared secret key, we have 
$$Curve(a,Curve(b,9))=Curve(b,Curve(a,9))$$ 
Then we further describe the details of shared secret computation and key generation. But before doing this, we first introduce two concept about the coordinate on the curve.

### Radix-$2^r$ Representation
A radix-$2^r$ representation represents an element $f$ in a $b$-bit prime field as $(f_0, f_1, \ldots, f_{\lceil b/r \rceil - 1}),$ s.t. 
$$f = \sum_{i=0}^{\lceil b/r \rceil - 1} f_i \cdot 2^{ir}$$
This is called the Radix-$2^r$ representation. Field arithmetic can then be carried out on limbs. Suppose there are 4 limbs and we want to do the simple addition of $a+b$, we can perform $a[i]+b[i]$, where $i=0,1,2,3$, and add corresponding limbs of the operands. Now we analyze two instances of radix-$2^r$ representation.
#### Radix-$2^{51}$ Representation
We represent an interger $f$ modulo $2^{255}-19$ as 
$$f_0 + 2^{51}f_1 + 2^{102}f_2 + 2^{153}f_3 + 2^{204}f_4$$.
The result of product of $f_0 + 2^{51}f_1 + 2^{102}f_2 + 2^{153}f_3 + 2^{204}f_4$ and $g_0 + 2^{51}g_1 + 2^{102}g_2 + 2^{153}g_3 + 2^{204}g_4$ is $h_0 + 2^{51}h_1 + 2^{102}h_2 + 2^{153}h_3 + 2^{204}h_4$ modulo $2^{255}-19$ where
$$h_0 = f_0g_0 + 19f_1g_4 + 19f_2g_3 + 19f_3g_2 + 19f_4g_1$$

$$h_1 = f_0g_1 + \quad f_1g_0 + 19f_2g_4 + 19f_3g_3 + 19f_4g_2$$

$$h_2 = f_0g_2 + \quad f_1g_1 + \quad f_2g_0 + 19f_3g_4 + 19f_4g_3$$

$$h_3 = f_0g_3 + \quad f_1g_2 + \quad f_2g_1 + \quad f_3g_0 + 19f_4g_4$$

$$h_4 = f_0g_4 + \quad f_1g_3 + \quad f_2g_2 + \quad f_3g_1 + \quad f_4g_0$$
This radix-$2^{51}$ is to fit the $64 \times 64 \rightarrow 128\text{-bit}.$ serial multiplier.

#### Radix-$2^{25.5}$ Representation
We represent an integer $f$ modulo $2^{255}-19$ as
$$f_0 + 2^{26}f_1 + 2^{51}f_2 + 2^{77}f_3 + 2^{102}f_4 + 2^{128}f_5 + 2^{153}f_6 + 2^{179}f_7 + 2^{204}f_8 + 2^{230}f_9$$
which means a 255-bits consists 10 limbs, five are 25 bits long and five are 26 bits long. Here we omit the result of $h=f \times g \quad mod 2^{255}-19$ since the equations are very similar to the previous radix-$2^{51}$ representation.

This representation is designed to fit the $32 \times 32 \rightarrow 64\text{-bit}$ vector multiplier, which can perform two ($128/64=2$) multiplications in one instruction. We can compute $h=f \times g$ and $h' = f'g'$ at the same time.

Then we will introduce an important part of how to delay the carry propagation. We perform a carry from $h_k$ to $h_{k+1}$ (the indice is in $\mathbb{F}_{10}$) in 3 steps:
- Perform logical right shift for the $64\text{-bit}$ word in $h_k$. The shifting amount is $26-(k mod 2)$
- Add the result from step one to $h_{k+1}$
- Mask the most significant $38+(kmod2)$ bits of $h_k$ i.e. set them to $0$. The head room of the unused space is for stroing carry bits, which is the key character of delaying carry propagation. 

Here we simply illustrate the carry propagation chain.
$$ h_0 \rightarrow h_1 \rightarrow h_2 \rightarrow h_3 \rightarrow h_4 \rightarrow h_5 \rightarrow h_6 \rightarrow h_7 \rightarrow h_8 \rightarrow h_9 \rightarrow h_0 \rightarrow h_1.$$

### Comparision with AVXECC
The operands used by AVX2 vector multiply instruction have to be stored in $256\text{-bit}$ register. $128 \times 128 \rightarrow 256$. so we can execute $32\text{-bit}$ or $64\text{-bit}$ elements in parallel. In AVXECC we execute $4 \times 1$ operation, which leaves maximum $32\text{-bit}$ length for each limbs. We leave the $3\text{-bit}$ head room to delay carry bits propagation, i.e. radix-$2^{29}$. 
 
### Affine Coordinates
Affine coordinates are the standard way to represent points in a two-dimensional space. A point $P$ defined on a elliptic curve $E$ defined by the equation $y^2=x^3 + ax +b$ in affine coordinates should satisfy this equation in the field.

More on extended affine coordinates. Beyond standard affine coordinate $(x,y)$, we use extended affine coordinates to deal with edge cases and ensure const time operations to counter side-channel attacks. For a Twisted Edwards curve
$$ax^2+y^2=1+dx^2y^2$$
the point $P$ on this curve can be represented as $(x,y,z,t)$, where 
$$x=X/Z,\quad y=Y/Z, \quad t=XY/Z$$

### Projective Coordinates
Projective coordinates is a alternative way to represent points on elliptic curve. In this form we can do calculation more efficiently. A point P in projective coordinates form is represented by three coordinates $(X:Y:Z)$. The relationship between projective coordinates $(X:Y:Z)$ and affine coordinates $(x,y)$ is given as follows:
$$x=X/Z,\quad y=Y/Z,\quad Z \neq 0$$ 
Why we use projective coordinates? In $\mathbb{F}_p$ with large $p$, divisions are a lot more expensive than multiplications. We use projective coordinates to delay inversions and divisions. 


### Shared-Secret Computation
The compuatation is x-coordinate only variable base scalar multiplication on Montgomery Curve, so we apply montgomery ladder, which is the most efficient algorithm.  For each bit of the scalar the ladder performs a differential addition and a doubling, which is call the **ladder step**. Note that the montgomery ladder use projective coordinates.



## A Lightweight EdDSA Signature Verification

In this section, we will go through the important part of the paper, which proposed a faster algorithem to compute two-scalar multiplication and implements it into EdDSA signature verfication. At the verification phase, we compute 
$$R=sB-hA $$ where $s$, $h$ are scalar, $B$ is the fixed point on the elliptic curve and $A$ is the secret key. This computation is a conbination of a fixed-base scalar multiplication and a variable-base multiplication, which can be done seperatly to gain higher efficience. We do the fixed-base compuation on Edwards25519 and the variable-base computaion on Curve25519. As we know these two curves are birational equevalant, so we can achieve the seperate scalar multiplication by point conversion.  
First we would like to introduce serveral related concepts to better understand this method.

### Non-Adjacent Form (NAF)

The NAF is a low-weight form to represent scalar. Each non-zero digit in the NAF is odd and is in the range of $[-2^{w-1},2^{w-1}]$, where $w$ is the width of the window and is a positive integer. We want the NAF to contain as many zero as possible, which is the key character of its efficiency. This method can reduce the number of required operation when performing multiplication. An integer $k$ can be represented in NAF as follows:
$$ k = \sum_{i=0}^{l} n_i 2^i$$ 
where $n_i$ is in $[-1,0,1]$, and no two consecutive $n_i$ is non-zero. 

### Joint-Sparse Form (JSF)
The JSF is a method to represent pairs of integers $(k,l)$ as a sum of power of $2$. The key character of JSF is that it is a sparse representation (meaning that it tries to minimize the non-zero terms in the sum) and compute the two scalar multiplication simultaneously. A pair of integers $(k,l)$ is represented by two sequence $(k_i,l_i)$ s.t. 
- each $k_i$ and $l_i$ is in the set $\{-1,0,1\}$. 
- $k_i$ and $l_i$ cannot be zero at the same time for any $i$.
- There are no two consecutive non-zero terms.

When computing $kP+lQ$, where $P$ and $Q$ are points on elliptic curve, the simultaneous processing of both scalars allows for the use of precomputed points $\pm P, \pm Q, \pm (P \pm Q)$.

### Simultaneous Double-Scalar Multiplications

This section introduces the method of computing $R=sB-hA$ using the JSF form. We skip this part since we prove that the separate scalar multiplication by expoliting the conversion between Twisted Edwards curve and Montgomery Cuvre is more efficient.

### Two Separate Scalar Multiplications

Now we can discuss the method to compute $R=sB-hA$ efficiently. We split the computation into two part i.e the fixed-base part $sB$ and the variable-base part $hA$. 

In this paper, it applys montgomery ladder to compute the fixed-base multiplication. For further improvement, we can use the look up table method in AVXECC. The computaion is originally done on TE curve, so we don't need to do the costly conversion if we use the look-up table method. And for the variable-base computation $hA$, we map it to montgomery curve and apply the montgomery ladder method. 

We can make the following improvements:
1. In this paper, it applys montgomery ladder to compute the fixed-base multiplication. We can use the look up table method in AVXECC.
2. We can extend the acceleration to key generation phase and signature generation phase. In key generation, there is one step we need to comupte the public key $sB$, where s is a scalar and B is a fixed point on the TE curve. So this is a fixed-base multiplication that can be optimized. In the signature generation phase, we need to compute $R=rB$, where r is a scalar and B is a fixed point. 
3. We can apply the SIMD style programming and parallel scalar multiplication to accelerate the fixed-base multiplication in key generation, signature generation and the two scalar multiplication in verification, which can be computed separately as one fixed-base multiplication and one varibale-base multiplication. In this way we can achieve high-throughput EdDSA signature implementaion. 

## Montgomery ladder and ladder step

This is a simple example of the ladder:

```python

def cswap(bit, R, S): # this is the constant time conditional swap
    dummy = bit * (R-S) # dummy = 0 or R-S
    R = R - dummy # R = R or R = S
    S = S + dummy # S = S or S = R
    return (R,S)

a = 44444 # this is the secret scalar
l = max # the maximum bit length, matching the order p
A = a.digits(2, padto = l) # pad a to the max length l
P0 = 0 # initial doubling 0 = 0P
P1 = P # difference P1 - P0 = P
for i in range (l-1, -1, -1): #from the most significant bit
    (P0,P1) = cswap(A[i], P0, P1)
    P1 = P0 + P1 #addition with fixed difference
    P0 = 2*P0 
    # here why we have this order? because we always 
    # want to overide the point we will not use first and the compute the doubling.
    (P0,P1) = cswap(A[i], P0, P1)
print(P0) # the results is stored in P0
# if the bit is zero, we double the smaller point the P0 and we add into the P1.
# if the bit is one, we double the bigger point and add into the lower one. 
```
If A[i] = 0, the new values are $P_0=2P_0,P_1=P_0+P_1$, if A[i] = 1, the new values are (after swaping back) $P_1=2P_1,P_0=P_0+P_1$. If the bit is zero, we double the smaller point the P0 and we add into the P1. If the bit is one, we double the bigger point and add into the lower one. Here either way the difference of $P_1-P_0=P$ is fixed. Addition is of points with known difference is called **differential addition**.

Here we give a simple numercial example of how the mon ladder works. Supppose the scalar is 9, we want to compute $9\times P$. We decode the scalar into binary code 1001, we set $P_0=0, P_1=P$. Start from the most significant bit.
- Iteration 1: A[3]=1, so we swap $P_0, P_1$, now $P_0=P, P_1=0$, we overide the point we will not use i.e. $P_1=P_1+P_0=P$, then we do the doubling $P_0=2 \times P_0=2P$. We swap them back, so after first iteration, we get $P_0=P, P_1=2P$.
- Iteration 2: A[2]=0, $P_1=P_1+P_0=3P$, $P_0=2\times P_0=2P$.
- Iteration 3: A[1]=0, $P_1=P_1+P_0=5P$, $P_0=2\times P_0=4P$.
- Iteration 4: A[0]=1, so we swap $P_0, P_1$. $P_1=P_1+P_0=9P$, $P_0=2\times P_0=10P$. After swap back, $P_0=9P, P_1=10P$.
- Return the result $P_0=9P$, which is exactly the result we want.

### Montgomery differential addition

Let $nP=(U_n:V_n:Z_n), mP=(U_p:V_p:Z_p)$ with known difference $(m-n)P=(U_{m-n}:V_{m-n}:Z_{m-n})$ on 
$$M_{A,B}:Bv^2=u^3+Au^2+u$$
We only use $X$ and $Z$.

- Addition: $n \neq m$
$$U_{m+n} = Z_{m-n} \left( (U_m - Z_m)(U_n + Z_n) + (U_m + Z_m)(U_n - Z_n) \right)^2$$
$$Z_{m+n} = U_{m-n} \left( (U_m - Z_m)(U_n + Z_n) - (U_m + Z_m)(U_n - Z_n) \right)^2$$
- Doubling: $n=m$
$$4U_nZ_n = (U_n + Z_n)^2 - (U_n - Z_n)^2$$
$$U_{2n} = (U_n + Z_n)^2(U_n - Z_n)^2$$
$$Z_{2n} = 4U_nZ_n \left( (U_n - Z_n)^2 + \left( \frac{A + 2}{4} \right)(4U_nZ_n) \right)$$

Differential addtion takse **4M** and **2S**, doubling takes **3M** and **2S**. In the ladder, $m-n=1$, choose $Z_{m-n}=1$ and $(A+2)/4$ small. Then cost per bit is **5M** and **4S**.

I will add the details of the ladder step soon...

## High-speed High-security signatures
We only focus on the fast software implementation of signatures, including: 
- fast signature verification, single or batched.
- fast signing.
- fast key generation

We would like to transfer such idea to fast EdDSA signature. First we will have some recap on EdDSA keys and signatures.

### EdDSA keys and signatures
An EdDSA secret key is a $b\text{-bit}$ string $k$. the hash (SHA-512) $H(k)=(h_0,h_1,...,h_{2b-1})$ determines an interger
$$a = 2^{b-2} + \sum_{3 \leq i \leq b-3} 2^{hi} \in \{2^{b-2}, 2^{b-2} + 8, \ldots, 2^{b-1} - 8\},$$
which in turn determines the multiple $A=aB$, and $A$ is the corresponding public key.
The signature of message $M$ under this secret key $k$ is defined as follows:
- Define $r=H(M)$, $R=rB$.
- Define $S=(r+H(R,A,M))\text{mod}l$.
The signature of $M$ under $k$ is then the $2b\text{-bit}$ string $(R,S)$, where $S$ is the $b\text{-bit}$ little-endian encoding of $S$.

Verification of the signature of $M$ under public key $A$ is defined as follows:
- Verifier computes $H(R,A,M)$ and then check the group equation $8SR=8R+8H(R,A,M)$

### Mitigating Bottleneck

For the signature generation part, the bottleneck is computing $R=rB$ for given $r$ and fixed point $B$. We write the scalar $r$ in the form of 
$$\sum_{i=0}^{63} 16^{i}r_i, \text{ where } r_i \in \{-8, -7, \ldots, 6, 7\}.$$
Finally we compute $rB$ as $\sum_i 16^{i}r_iB.$
