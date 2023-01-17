# Nathaniel Keninger
# Creation Date 10/23/2022

# Tg Calculations for Liquid Fragility Data, based on work done by Mauro et. al.

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
import csv
import warnings




#########################################################################
#                       Lithium Borate Structure                        #
#########################################################################

i = np.arange(0, 3, .001)

N4 = []
B3 = []
T3 = []
T4 = []
D3 = []
D4 = []
L4 = []
M3 = []
P3 = []
O3 = []


for R in i:
    if R < .25:
        n4 = R
        b3 = 1 - 4*R
        d3 = 0
        d4 = 0
        t3 = 3*R
        t4 = t3/3
        l4 = 0
        m3 = 0
        p3 = 0
        o3 = 0       
    if .25 <= R < .4:
        n4 = R
        b3 = 0
        d3 = 2*R - .5
        d4 = d3
        t3 = 3*(R - d3)
        t4 = t3/3
        l4 = 0
        m3 = 0
        p3 = 0
        o3 = 0
    if .4 <= R < .7:
        n4 = 1/3 + 1/6*R
        b3 = 0
        d3 = .5*(1 - R)
        d4 = d3
        t3 = 1 - R - d3
        t4 = t3/3
        l4 = .5*R - (2/3) * d3
        m3 = 5/6*R - 1/3
        p3 = 0
        o3 = 0
    if .7 <= R < 1:
        n4 = 5/8 - 1/4*R
        b3 = 0
        d3 = .5*(1-R)
        d4 = d3
        t3 = 1 - R - d3
        t4 = t3/3 
        l4 = 7/24 + 1/12*R - 2/3*d3
        m3 = 5/4*R - 5/8
        p3 = 0
        o3 = 0
    if 1 <= R < 2:
        n4 = 5/8 - 1/4*R
        o3 = ((1/6)/.86)*R - ((1/6)/.86)
        b3 = 0
        d3 = 0
        d4 = d3
        t3 = 0
        t4 = t3/3 
        l4 = n4
        m3 = 11/8 - .75*R + o3
        p3 = R - 1 - 2*o3
        R_2 = R
        p3_2 = p3
        m3_2 = m3
        o3_2 = o3
    if 2 <= R < 2.5:
        n4 = 5/8 - 1/4*R
        b3 = 0
        d3 = 0
        d4 = 0
        t3 = 0
        t4 = 0
        l4 = n4
        p3 = ((.5 - p3_2)/.5)*(R - 2) + p3_2#-.2248*R + 1.062
        o3 = ((.5 - o3_2)/.5)*(R - 2) + o3_2
        m3 = 1 - o3 - l4 - p3
        R_2_5 = R
    if 2.5 <= R < 3:
        n4 = 0
        b3 = 0
        d3 = 0
        d4 = 0
        t3 = 0
        t4 = 0
        l4 = 0
        m3 = 0
        p3 = 3 - R
        o3 = R - 2
        
    N4.append(n4)    
    B3.append(b3)
    D3.append(d3)
    D4.append(d4)
    T3.append(t3)
    T4.append(t4)
    L4.append(l4)
    M3.append(m3)
    P3.append(p3)
    O3.append(o3)

units = [B3, D3, D4, T3, T4, L4, M3, P3, O3]
unitNames = ['B3', 'D3', 'D4', 'T3', 'T4', 'L4', 'M3', 'P3', 'O3']

def testPlot(B3):
    plt.plot(i, B3)

    plt.show()


def liPlot():
    for unit in range(0,len(units)):
        if units[unit] != N4:
            plt.plot(i, units[unit], label = str(unitNames[unit]))
    

    plt.legend()
    

    plt.show()

#########################################################################
#                   Lithium Borate Network Formers                      #
#########################################################################


Ω = []
B4c = []
D4c = []
B3c = []
OB = []
LiNB = []
LiL = []

ΩNoD = []
B4cNoD = []
D4cNoD = []
B3cNoD = []
OBNoD = []
LiNBNoD = []

netForm = [B4c, B3c, D4c, OB, LiNB, LiL]
netFormName = ['B4*', 'B3*', 'D4*', 'OB', 'LiNB', 'LiLoose'] 

for R in range(0,len(i)):
    looseM = 2 * i[R] - L4[R] - M3[R] - P3[R] - 1.5*O3[R]
    ω = L4[R] + T4[R] + D4[R] + D3[R] + B3[R] + T3[R] + M3[R] + (2*(T4[R] + D4[R] + L4[R])) + 1.5*(B3[R] + D3[R] + T3[R]) + M3[R] + M3[R] + L4[R] + P3[R] + 1.5*O3[R] + looseM
    
    b4c = (L4[R] + T4[R])/ω
    b3c = (B3[R] + T3[R] + D3[R] + M3[R])/ω
    d4c = (D4[R])/ω
    ob = (2*(T4[R] + D4[R] + L4[R]) + 1.5*(B3[R] + D3[R] + T3[R]) + M3[R])/ω
    linb = (L4[R] + M3[R] + P3[R] + 1.5*O3[R])/ω
    liloose = looseM/ω

    ωNoD = ω - looseM

    b4cNoD = (L4[R] + T4[R])/ωNoD
    b3cNoD = (B3[R] + T3[R] + D3[R] + M3[R])/ωNoD
    d4cNoD = (D4[R])/ωNoD
    obNoD = (2*(T4[R] + D4[R] + L4[R]) + 1.5*(B3[R] + D3[R] + T3[R]) + M3[R])/ωNoD
    linbNoD = (L4[R] + M3[R] + P3[R] + 1.5*O3[R])/ωNoD
    

    
    Ω.append(ω)
    B4c.append(b4c)
    B3c.append(b3c)
    D4c.append(d4c)
    OB.append(ob)
    LiNB.append(linb)
    LiL.append(liloose)

    ΩNoD.append(ωNoD)
    B4cNoD.append(b4cNoD)
    B3cNoD.append(b3cNoD)
    D4cNoD.append(d4cNoD)
    OBNoD.append(obNoD)
    LiNBNoD.append(linbNoD)
    

def netFormPlot():
    for unit in range(0, len(netForm)):
        plt.plot(i, netForm[unit], label = netFormName[unit])

    plt.legend()
    plt.show()



#########################################################################
#                         Constraint Counting                           #
#########################################################################



α = []
β = []
μ = []
γ = []
δ = []

nTg = []
nT = []

αNoD = []
βNoD = []
μNoD = []
γNoD = []

nTgNoD = []
nTNoD = []

constraints = [α, β, μ, γ, δ]
constraintName = ['α', 'β', 'μ', 'γ', 'δ']

for R in range(0, len(i)):
    a = OB[R] * 2
    b = 5*B4c[R] + 4.5*D4c[R] + 3*B3c[R]
    m = 2*LiNB[R]
    g = OB[R] #+ LiNB[R]
    d = 2.5*LiL[R]

    aNoD = OBNoD[R] * 2
    bNoD = 5*B4cNoD[R] + 4.5*D4cNoD[R] + 3*B3cNoD[R]
    mNoD = 2*LiNBNoD[R]
    gNoD = OBNoD[R] #+ LiNB[R]
        
    ntg = a + b + m  + d
    nt = a + b + m + g + d

    ntgNoD = aNoD + bNoD + mNoD
    ntNoD = aNoD + bNoD + mNoD + gNoD 

    
    α.append(a)
    β.append(b)
    μ.append(m)
    γ.append(g)
    δ.append(d)

    nTg.append(ntg)
    nT.append(nt)


    αNoD.append(aNoD)
    βNoD.append(bNoD)
    μNoD.append(mNoD)
    γNoD.append(gNoD)

    nTgNoD.append(ntgNoD)
    nTNoD.append(ntNoD)




    
def partConstraintPlot():
    for c in range(0, len(constraints)):
        plt.plot(i, constraints[c], label = constraintName[c])
    plt.legend()
    plt.show()

def netConstraintPlotTg():
    plt.plot(i, nTg)
    plt.show()

def netConstraintPlot():
    plt.plot(i, nT)
    plt.show()
    
#########################################################################
#                            Tg Calculation                             #
#########################################################################

Tg_r = 533 #K , Reference Tg for Lithium Borate Glass System
n_r = 2.4 # Reference Composition for Borate Glass System

rDataAffatig = [0,.05,.1,.15,.2,.2,.3,.3,.4,.5,.6,.7,.8,.9,1,1.2,1.5]
TgDataAffatig = [533,570,616,668,708,706,757,762,769,773,756,741,730,714,694,651,602]

rDataWhiteBook = [1.5,1.9,2,2.2,2.4,2.6,2.7,2.75,2.75,2.8,3]
TgDataWhiteBook = [604,566,554,547,542,543,544,548,544,546,568]

xDataKodama = [.00,.02,.04,.06,.08,.1,.12,.14,.16,.18,.2,.22,.24,.26,.28]
RDataKodama = []
for j in range(0, len(xDataKodama)):
    rdatakodama = xDataKodama[j]/(1-xDataKodama[j])
    RDataKodama.append(rdatakodama)
TgDataKodama = [518,528,543,558,598,618,638,673,693,723,748,753,758,763,763]


Tg = []

TgNoD = []
for R in range(0, len(i)):
    tg = (3 - n_r)/(3-nTg[R]) * Tg_r
    tgNoD = (3 - n_r)/(3-nTgNoD[R]) * Tg_r

    Tg.append(tg)

    TgNoD.append(tgNoD)

TgMauro = []

for R in i:
    x = R/(1+R)
    if R <= .5:
        tgmauro = (1/5)*((5 - 4*x)/(1 - 2*x))*543
    elif R <= 1:
        tgmauro = (1/5)*((31 - 38*x)/(8*x - 1))*543
    else:
        tgmmauro = 0
    if tgmauro > 760:
        tgmauro = 760
        
    TgMauro.append(tgmauro)


def TgPlot():
    plt.plot(i, Tg, label = "Tg Model")
    plt.plot(rDataAffatig, TgDataAffatig, 'o',c = "red", label = "Affatigato et. al.")
    plt.plot(RDataKodama, TgDataKodama, '^',c = "green", label = "Kodama and Kojima")
    plt.plot(rDataWhiteBook, TgDataWhiteBook, 'o', c = "purple", label = "White Book")

    plt.xlim(0, 3)
    plt.ylim(450, 850)
    plt.title("Glass Transition Temperature for Lithium Borate Glass")
    plt.xlabel("R")
    plt.ylabel("Tg(K)")
    
    plt.legend()
    plt.show()
    

def TgPlotNoD():
    plt.plot(i, TgNoD, label = "Takeda Model")
    plt.plot(rDataAffatig, TgDataAffatig, 'o',c = "red", label = "Affatigato et. al.")
    plt.plot(RDataKodama, TgDataKodama, '^',c = "green", label = "Kodama and Kojima")
    plt.plot(rDataWhiteBook, TgDataWhiteBook, 'o', c = "purple", label = "White Book")
    
    plt.xlim(0, 3)
    plt.ylim(450, 850)
    plt.title("Tg for Lithium Borate Glass")
    plt.xlabel("R")
    plt.ylabel("Tg(K)")

    plt.legend()
    plt.show()


def TgPlotVsNoδ():
    plt.plot(i, Tg, label = "δ-constraint Model")
    plt.plot(i, TgNoD, label = "Takeda Model")
    plt.plot(i, TgMauro, label = "Mauro-Gupta Model")
    plt.plot(rDataAffatig, TgDataAffatig, 'o',c = "red", label = "Affatigato et. al.")
    plt.plot(RDataKodama, TgDataKodama, '^',c = "green", label = "Kodama and Kojima")
    plt.plot(rDataWhiteBook, TgDataWhiteBook, 'o', c = "purple", label = "White Book")

    plt.xlim(0, 3)
    plt.ylim(450, 850)
    plt.title("Tg for Lithium Borate Glass with and without δ-constraint effects")
    plt.xlabel("R")
    plt.ylabel("Tg(K)")

    plt.legend()
    plt.show()



#########################################################################
#                         Fragility Calculation                         #
#########################################################################

k = sc.k # 1.380649e-23 m^2 kg s^-1 K^-1 (Boltzmann cst)

m_0 = 14.97 # fragility of a strong glass

# BORATE GLASS SYSTEM CONSTRAINT ONSET TEMPERATURES

T_α = 921 # from previous lit (collin's YM paper)
T_β = 715 # needs to be significantly above Tg such that discrete model applies for all Tg
T_γ = 393 # from previous lit(collin's YM paper)
T_μ = 900 # above Τ_β and below T_α

T_δ = 850 # likely near or below T_μ my guess is slightly below



vt = 750 # Escape attemps, fitting parameter fit to first experimental point

def ΔF_c(T_c):   # Activation free energy for constraint 'c'
    ΔF_c = -k * T_c * np.log(1 - 2**(-1/(vt)))
    return ΔF_c

def q(T_c,T):    # Temperature dependence of given constraint
    q = (1 - np.exp(-ΔF_c(T_c)/(k*T)))**vt
    return q

def f(R, T):
    f = 3 - (α[R] * q(T_α, T)) - (β[R] * q(T_β, T)) - (μ[R] * q(T_μ, T)) - (γ[R] * q(T_γ, T)) - (δ[R] * q(T_δ, T))
    return f

def fNoD(R, T):
    f = 3 - (αNoD[R] * q(T_α, T)) - (βNoD[R] * q(T_β, T)) - (μNoD[R] * q(T_μ, T)) - (γNoD[R] * q(T_γ, T))
    return f

m = []
mδ = []

d = 1

for R in range(0, (len(i))):
    if R < len(i) - d:
        T = TgNoD[R]
        T2 = TgNoD[R + d]
        F = fNoD(R,T)
        F2 = fNoD(R, T2)

        warnings.filterwarnings("ignore")
        
        df = (np.log(F2) - np.log(F))
        dT = (np.log(T2) - np.log(T))
        dfdT = df/dT


        M = m_0*(1 + dfdT)

        T = Tg[R]
        T2 = Tg[R + d]
        F = f(R,T)
        F2 = f(R, T2)


        df = (np.log(F2) - np.log(F))
        dT = (np.log(T2) - np.log(T))
        dfdT = df/dT

        Mδ = m_0*(1 + dfdT)
    else:
        M = 0
        Mδ = 0
    m.append(M)
    mδ.append(Mδ)



# Data Points
# Chyrssikos et al, Table II of Mauro Fragility Paper
xChrys = [0.00,0.02,0.06,0.08,0.1,.12,.16,.18,.19,.20,.21,.25]

RChrys = []

for j in range(0, len(xChrys)):
    rchrys = xChrys[j]/(1-xChrys[j])

    RChrys.append(rchrys)

mChrys = [33,35,37,36,37,37,44,47,51,61,55,61]

ix = []

for R in i:
    
        
    ixval = R/(1+R)

    ix.append(ixval)


def fragilityPlot():
    plt.plot(i, m, label = "TCT Model")
    plt.plot(RChrys, mChrys, 'o', label = "Experimental Data")
    
    plt.xlim(0, 1)
    plt.ylim(0, 70)

    plt.title("Fragility of Lithium Borate Glass")
    plt.xlabel("R")
    plt.ylabel("Fragility")
    plt.legend()

    

    plt.show()
    
def fragilityPlotX():
    plt.plot(ix, m, label = "TCT Model")
    plt.plot(xChrys, mChrys, 'o', label = "Experimental Data")
    
    plt.xlim(0, .5)
    plt.ylim(0, 70)

    plt.title("Fragility of Lithium Borate Glass")
    plt.xlabel("x")
    plt.ylabel("Fragility")
    plt.legend()

    plt.show()

def fragilityPlotδ():
    plt.plot(i, m, label = "TCT Model")
    plt.plot(i, mδ, 'green', label = "TCT Model with δ-constraint")
    plt.plot(RChrys, mChrys, 'o', label = "Experimental Data")
    
    plt.xlim(0, 1)
    plt.ylim(0, 70)

    plt.title("Fragility of Lithium Borate Glass")
    plt.xlabel("R")
    plt.ylabel("Fragility")
    plt.legend()

    

    plt.show()

def fPlotTest():
    plt.plot(o, p)
    plt.plot(o,p)
    plt.show()
    

def exportFragility():

    head = ['x', 'R', 'm', 'md']
    data = []
    for c in range(0,len(i)):
        #line = {'x': ix[c], 'R': i[c], 'm': m[c], 'mδ': mδ[c]}
        line = [ix[c], i[c], m[c], mδ[c]]
        c+= 1
        data.append(line)

    with open("FragilityData.csv", 'w', newline = '') as file:
    
        #writer = csv.DictWriter(file, fieldnames = head)
        writer = csv.writer(file)
        
        writer.writerow(head)
        writer.writerows(data)

    with open("FragilityExpData.csv", 'w', newline = '') as file:
        writer = csv.writer(file)
        head = ['x', 'R', 'Chyrssikos et al']
        expData = []
        for a in range(0,len(xChrys)):
            line = [xChrys[a], RChrys[a], mChrys[a]]
            expData.append(line)
            
        writer = csv.writer(file)
        writer.writerow(head)
        writer.writerows(expData)
    print(data)

#########################################################################
#                      Young's Modulus Calculation                      #
#########################################################################

#E = dE/dΔFc[Fc(x) - Fc(x)]

#Fc(x) = ρ(x)N_a/M(x) * Σ_(l=constraints) ΔF_l * n_l(x) * q_l(T)

#dE/dΔFc and F`c(x) are fitting parameters

FC = 1900
dEdΔFc = .033
ρ = 2.2


M = []
M_B2O3 = 2*10.811 + 3*15.999
M_Li2O = 2*6.9410 + 15.999
for R in i:
    mm = R * M_Li2O + M_B2O3
    M.append(mm)



def F_c(R):
    Fc = (ρ*sc.N_A/M[R]) * (ΔF_c(T_α)*α[R]*q(T_α,T) + ΔF_c(T_β)*β[R]*q(T_β,T) + ΔF_c(T_γ)*γ[R]*q(T_γ,T) + ΔF_c(T_δ)*δ[R]*q(T_δ,T) + ΔF_c(T_μ)*μ[R]*q(T_μ, T))
 

def menu():
    print("="*50)
    print("="*50)
    choosing = True
    while choosing:
        try:
            a = int(input("Please select an operation:\n1. Lithium Borate Structure Plot\n2. Tg Plot\n3. Tg Plot Compared to Other Models\n4. Fragility Plot\n5. Fragility Plot with δ-constraint effects (experimental)\n6. Young's Modulus Plot\n0. Exit\n"))
        
            if a in range(1,7):           
                if a == 1:
                    liPlot()
                if a == 2:
                    TgPlot()
                if a == 3:
                    TgPlotVsNoδ()
                if a == 4:
                    fragilityPlot()
                if a == 5:
                    fragilityPlotδ()
                if a == 6:
                    print("Work In Progress, Try again later.")
                choosing = False
                menu()
            elif a == 0:
                choosing = False
            else:
                print("Please choose a valid number")
        except ValueError:
            print("Please choose a valid number")
    

menu()
