import time
import quantum

print("Provide the number below.")
print("1. Only RHF calculation.")
print("2. RHF + MP2")
print("3. RHF + MP2 + MP3")
print("4. RHF + MP2 + MP3 + LCCD")
print("5. RHF + MP2 + MP3 + CID")
print("6. RHF + MP2 + MP3 + LCCD + CID + finite-field")

control = input("Your choice: ")
control = int(control)

start_time = time.time()

if (control == 1): #RHF only
    start_hf_time = time.time()
    nirrep, n, ej, noc, nu, n2el, natms, coord, nch = quantum.driver()
    s12, ovrl, ekin, epot, dipole_raw = quantum.prefock(n, n2el) #Loading integrals from prefock procedure
    ca, escf, fap, pa = quantum.fock(ej, n, ekin, epot, noc, s12) #RHF calculations
    end_hf_time = time.time()
    print("RHF procedure time:", round(end_hf_time - start_hf_time, 2), "seconds")
    dipole_moment, total_dipole_moment, dipnu, dipel = quantum.dipole(natms, coord, nch, n, dipole_raw, ca, noc, pa) #Dipole moment calculations
    quadrupole_moment = quantum.quadrupole(natms, nch, coord, n, dipole_raw, pa) #Quadruple moment calculations
    end_time = time.time()
    print("Total calculations time: ", round(end_time - start_time, 2), "seconds")

elif (control == 2): #RHF + MP2
    start_hf_time = time.time()
    nirrep, n, ej, noc, nu, n2el, natms, coord, nch = quantum.driver()
    s12, ovrl, ekin, epot, dipole_raw = quantum.prefock(n, n2el) #Loading integrals from prefock procedure
    ca, escf, fap, pa = quantum.fock(ej, n, ekin, epot, noc, s12) #RHF calculations
    end_hf_time = time.time()
    print("RHF procedure time:", round(end_hf_time - start_hf_time, 2), "seconds")
    dipole_moment, total_dipole_moment, dipnu, dipel = quantum.dipole(natms, coord, nch, n, dipole_raw, ca, noc, pa) #Dipole moment calculations
    quadrupole_moment = quantum.quadrupole(natms, nch, coord, n, dipole_raw, pa) #Quadruple moment calculations
    print("Integrals transformation in progress...")
    start_trans_time = time.time()
    quantum.trans4(n, ca) #Transformation
    end_trans_time = time.time()
    print("Integrals transformation time:", round(end_trans_time - start_trans_time, 2), "seconds")
    vh, vp, vl, vr = quantum.sortab(n, noc, nu) #Sorting 
    print("MP2 calculations in progress...")
    c2, den, mp2 = quantum.mp2(noc, nu, fap, vr, escf) #MP2
    end_time = time.time()
    print("Total calculations time: ", round(end_time - start_time, 2), "seconds")
    
elif (control == 3): #HF + MP2 + MP3
    start_hf_time = time.time()
    nirrep, n, ej, noc, nu, n2el, natms, coord, nch = quantum.driver()
    s12, ovrl, ekin, epot, dipole_raw = quantum.prefock(n, n2el) #Loading integrals from prefock procedure
    ca, escf, fap, pa = quantum.fock(ej, n, ekin, epot, noc, s12) #RHF calculations
    end_hf_time = time.time()
    print("RHF procedure time:", round(end_hf_time - start_hf_time, 2), "seconds")
    dipole_moment, total_dipole_moment, dipnu, dipel = quantum.dipole(natms, coord, nch, n, dipole_raw, ca, noc, pa) #Dipole moment calculations
    quadrupole_moment = quantum.quadrupole(natms, nch, coord, n, dipole_raw, pa) #Quadruple moment calculations
    print("Integrals transformation in progress...")
    start_trans_time = time.time()
    quantum.trans4(n, ca) #Transformation
    end_trans_time = time.time()
    print("Integrals transformation time:", round(end_trans_time - start_trans_time, 2), "seconds")
    vh, vp, vl, vr = quantum.sortab(n, noc, nu) #Sorting 
    print("MP2 calculations in progress...")
    c2, den, mp2 = quantum.mp2(noc, nu, fap, vr, escf) #MP2
    print("MP3 calculations in progress...")
    mp3 = quantum.mp3(noc, nu, c2, vr, vl, vh, vp, escf, mp2) #MP3
    end_time = time.time()
    print("Total calculations time: ", round(end_time - start_time, 2), "seconds")

elif (control == 4): #HF + MP2 + MP3 + LCCD
    start_hf_time = time.time()
    nirrep, n, ej, noc, nu, n2el, natms, coord, nch = quantum.driver()
    s12, ovrl, ekin, epot, dipole_raw = quantum.prefock(n, n2el) #Loading integrals from prefock procedure
    ca, escf, fap, pa = quantum.fock(ej, n, ekin, epot, noc, s12) #RHF calculations
    end_hf_time = time.time()
    print("RHF procedure time:", round(end_hf_time - start_hf_time, 2), "seconds")
    dipole_moment, total_dipole_moment, dipnu, dipel = quantum.dipole(natms, coord, nch, n, dipole_raw, ca, noc, pa) #Dipole moment calculations
    quadrupole_moment = quantum.quadrupole(natms, nch, coord, n, dipole_raw, pa) #Quadruple moment calculations
    print("Integrals transformation in progress...")
    start_trans_time = time.time()
    quantum.trans4(n, ca) #Transformation
    end_trans_time = time.time()
    print("Integrals transformation time:", round(end_trans_time - start_trans_time, 2), "seconds")
    vh, vp, vl, vr = quantum.sortab(n, noc, nu) #Sorting 
    print("MP2 calculations in progress...")
    c2, den, mp2 = quantum.mp2(noc, nu, fap, vr, escf) #MP2
    print("MP3 calculations in progress...")
    mp3 = quantum.mp3(noc, nu, c2, vr, vl, vh, vp, escf, mp2) #MP3
    print("LCCD calculations in progress...")
    start_lccd_time = time.time()
    lccd, ecorr_lccd = quantum.lccd(c2, noc, nu, mp2, vr, vl, vh, vp, escf, den)
    end_lccd_time = time.time()
    print("LCCD calculations time: ", round(end_lccd_time - start_lccd_time, 2), "seconds")
    end_time = time.time()
    print("Total calculations time: ", round(end_time - start_time, 2), "seconds")

elif (control == 5): #HF + MP2 + MP3 + CID
    start_hf_time = time.time()
    nirrep, n, ej, noc, nu, n2el, natms, coord, nch = quantum.driver()
    s12, ovrl, ekin, epot, dipole_raw = quantum.prefock(n, n2el) #Loading integrals from prefock procedure
    ca, escf, fap, pa = quantum.fock(ej, n, ekin, epot, noc, s12) #RHF calculations
    end_hf_time = time.time()
    print("RHF procedure time:", round(end_hf_time - start_hf_time, 2), "seconds")
    dipole_moment, total_dipole_moment, dipnu, dipel = quantum.dipole(natms, coord, nch, n, dipole_raw, ca, noc, pa) #Dipole moment calculations
    quadrupole_moment = quantum.quadrupole(natms, nch, coord, n, dipole_raw, pa) #Quadruple moment calculations
    print("Integrals transformation in progress...")
    start_trans_time = time.time()
    quantum.trans4(n, ca) #Transformation
    end_trans_time = time.time()
    print("Integrals transformation time:", round(end_trans_time - start_trans_time, 2), "seconds")
    vh, vp, vl, vr = quantum.sortab(n, noc, nu) #Sorting 
    print("MP2 calculations in progress...")
    c2, den, mp2 = quantum.mp2(noc, nu, fap, vr, escf) #MP2
    print("MP3 calculations in progress...")
    mp3 = quantum.mp3(noc, nu, c2, vr, vl, vh, vp, escf, mp2) #MP3
    print("CID calculations in progress...")
    start_cid_time = time.time()
    cid, ecorr_cid = quantum.cid(c2, noc, nu, mp2, vr, vl, vh, vp, escf, den)
    end_cid_time = time.time()
    print("CID calculations time: ", round(end_cid_time - start_cid_time, 2), "seconds")
    end_time = time.time()
    print("Total calculations time: ", round(end_time - start_time, 2), "seconds")

elif (control == 6): #HF + MP2 + MP3 + LCCD + CID + finite-field properties (all)
    start_hf_time = time.time()
    nirrep, n, ej, noc, nu, n2el, natms, coord, nch = quantum.driver()
    s12, ovrl, ekin, epot, dipole_raw = quantum.prefock(n, n2el) #Loading integrals from prefock procedure
    ca, escf, fap, pa = quantum.fock(ej, n, ekin, epot, noc, s12) #RHF calculations
    end_hf_time = time.time()
    print("RHF procedure time:", round(end_hf_time - start_hf_time, 2), "seconds")
    dipole_moment, total_dipole_moment, dipnu, dipel = quantum.dipole(natms, coord, nch, n, dipole_raw, ca, noc, pa) #Dipole moment calculations
    quadrupole_moment = quantum.quadrupole(natms, nch, coord, n, dipole_raw, pa) #Quadruple moment calculations
    print("Integrals transformation in progress...")
    start_trans_time = time.time()
    quantum.trans4(n, ca) #Transformation
    end_trans_time = time.time()
    print("Integrals transformation time:", round(end_trans_time - start_trans_time, 2), "seconds")
    vh, vp, vl, vr = quantum.sortab(n, noc, nu) #Sorting 
    print("MP2 calculations in progress...")
    c2, den, mp2 = quantum.mp2(noc, nu, fap, vr, escf) #MP2
    print("MP3 calculations in progress...")
    mp3 = quantum.mp3(noc, nu, c2, vr, vl, vh, vp, escf, mp2) #MP3
    print("LCCD calculations in progress...")
    start_lccd_time = time.time()
    lccd, ecorr_lccd = quantum.lccd(c2, noc, nu, mp2, vr, vl, vh, vp, escf, den)
    end_lccd_time = time.time()
    print("LCCD calculations time: ", round(end_lccd_time - start_lccd_time, 2), "seconds")
    print("CID calculations in progress...")
    start_cid_time = time.time()
    cid, ecorr_cid = quantum.cid(c2, noc, nu, mp2, vr, vl, vh, vp, escf, den)
    end_cid_time = time.time()
    print("CID calculations time: ", round(end_cid_time - start_cid_time, 2), "seconds")
    start_ff_time = time.time()
    print("Finite-field calculations of electric properties in progress...")
    dipole_moment_ff_hf, total_dipole_moment_ff_hf, alpha_ff_hf, alpha_ave_ff_hf, dipole_moment_ff_mp2, total_dipole_moment_ff_mp2, alpha_ff_mp2, alpha_ave_ff_mp2, dipole_moment_ff_mp3, total_dipole_moment_ff_mp3, alpha_ff_mp3, alpha_ave_ff_mp3 = quantum.finite_field(n, dipole_raw, ekin, epot, ej, dipnu, noc, s12, nu, vr, vl, vh, vp, escf, mp2, mp3) #Własności elektryczne metodą skończonego pola
    end_ff_time = time.time()
    print("FF calculations time: ", round(end_ff_time - start_ff_time, 2), "seconds")
    end_time = time.time()
    print("Total calculations time: ", round(end_time - start_time, 2), "seconds")

else:
    print("You provided the wrong number! Try again.")


