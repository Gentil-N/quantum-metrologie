using QuantumOptics
using PyPlot
using LinearAlgebra

println("import done")

# Parameters
hbar = 1.0
N_photons = 2 #photon cutoff
hilbert_dim = N_photons + 1
N_atoms = 1 #number of atoms
cf = 1.0
af = 1.0
coupling = 1 #820.0
kappa = 1 #800.0 * 1000
repump = 1 #7.5 * 1000


# Bases
b_fock = FockBasis(hilbert_dim)
b_spin = SpinBasis(N_atoms//2)

# Fundamental operators
op_a = destroy(b_fock)
op_ad = create(b_fock)
op_n = number(b_fock)

op_sm = sigmam(b_spin)
op_sp = sigmap(b_spin)
op_sz = sigmaz(b_spin)


# Jaynes-Cummings-Hamiltonian
op_h_base = hbar * (cf - af) * op_n ⊗ sparse(one(b_spin))
op_hint = hbar * coupling / sqrt(N_atoms) * (op_ad ⊗ op_sm + op_a ⊗ op_sp)
op_h = op_h_base + op_hint


# Initial state
psi_init = fockstate(b_fock, 1) ⊗ spindown(b_spin)
#psi_init = coherentstate(b_fock, 1) ⊗ spindown(b_spin)

# Time interval
T_end= 5#0.1
dt= 0.05 #0.0001
time_range = [0:dt:T_end;]

# Collapse operators and decay rates
op_colapse = [sparse(one(b_fock)) ⊗ op_sp, op_a ⊗ sparse(one(b_spin))]
rates = [repump, kappa]

# Time evolution according to a master equation
print("computing...")
tout, result = timeevolution.master(time_range, psi_init, op_h, op_colapse; rates = rates)
println("done")

function get_u_from_h(t, h)
    return exp((-1im) * h * t / hbar)
end

print("g²...")
my_g2 = []
i = 1

function compute_g2(g2, index, time_range, op_a, op_ad, op_h, b_fock, b_spin)
    for t in time_range
        println(t)
        if t == 0.0
            op_u = sparse(one(b_fock)) ⊗ sparse(one(b_spin))
        else
            op_u = get_u_from_h(t, op_h)
        end
        op_ad_evol = dagger(op_u) * (op_ad ⊗ sparse(one(b_spin))) * op_u
        op_a_evol = dagger(op_u) * (op_a ⊗ sparse(one(b_spin))) * op_u
        norm_factor = (tr(result[index] * (op_ad ⊗ sparse(one(b_spin))) * (op_a ⊗ sparse(one(b_spin)))))^2
        if norm_factor == 0.0
            push!(g2, 0.0)
            index += 1
        end
        num = tr(result[index] * (op_ad ⊗ sparse(one(b_spin))) * op_ad_evol * op_a_evol * (op_a ⊗ sparse(one(b_spin)))) / norm_factor
        push!(g2, num)
        index += 1
    end
end

compute_g2(my_g2, i, time_range, op_a, op_ad, op_h, b_fock, b_spin)
println("done")


state_g0 = fockstate(b_fock, 0) ⊗ spindown(b_spin)
state_e0 = fockstate(b_fock, 0) ⊗ spinup(b_spin)
state_e1 = fockstate(b_fock, 1) ⊗ spinup(b_spin)
res_g0 = real(expect(state_g0 ⊗ dagger(state_g0), result))
res_g1 = real(expect(psi_init ⊗ dagger(psi_init), result))
res_e0 = real(expect(state_e0 ⊗ dagger(state_e0), result))
res_e1 = real(expect(state_e1 ⊗ dagger(state_e1), result))


figure(num=0)
plot(time_range, res_g0)
plot(time_range, res_g1)
plot(time_range, res_e0)
plot(time_range, res_e1)
xlabel(L"time")
ylabel(L"probabilities")

figure(num=1)
plot(time_range, real(my_g2))
xlabel(L"time")
ylabel(L"g²")
show()