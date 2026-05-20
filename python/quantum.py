import numpy as np
import sys

np.set_printoptions(legacy='1.25')


def driver():
    try:
        jobpar = open("jobpar", "r")
    except FileNotFoundError:
        print("jobpar file not found! Are you sure you provided the correct basis set and it is available in GENBAS?")
        sys.exit("Aborting execution!")
    nirrep = np.fromfile(jobpar, dtype="int64", count=1)  # number of irreducible representations
    nirrep = np.int64(nirrep[0])
    n = np.fromfile(jobpar, dtype="int64", count=60)  # number of basis functions
    n = np.int64(n[-1])
    ej = np.fromfile(jobpar, dtype="float64", count=60)  # nuclei repulsion energy
    ej = np.float64(ej[-1])
    noc = np.fromfile(jobpar, dtype="int64", count=60)  # number of electrons
    noc = np.int64(noc[-1])
    natms = np.fromfile(jobpar, dtype="int64", count=60)  # number of atoms
    natms = np.int64(natms[-1])
    coord_raw = np.fromfile(jobpar, dtype="float64", count=(59 + natms * 3))
    nch_raw = np.fromfile(jobpar, dtype="float64", count=(59 + natms * 3))
    n2el_raw = np.fromfile(jobpar, dtype="int64", count=(-1))  # number of non-zero two-electron integrals
    n2el = np.int64(0)
    for i in range(0, len(n2el_raw)):
        if n2el_raw[i] != 0:
            n2el = np.int64(n2el_raw[i])
            break

    k = 59
    coord = np.zeros((natms, 3))
    for j in range(0, natms):
        for i in range(0, 3):
            coord[j, i] = coord_raw[k]  # atoms' coordinates
            k = k + 1
    coord = np.transpose(coord)

    nch = []
    for i in range(0, len(nch_raw)):
        if nch_raw[i] != 0:
            nch.append(nch_raw[i])  # atomic charges
        else:
            continue

    print("Number of irreducible representations:", nirrep)
    print("Basis functions:", n)
    print("Nuclear repulsion energy:", ej)
    print("Number of electrons:", noc * 2)
    print("Number of non-zero two-electron integrals:", n2el)
    print("Number of atoms:", natms)
    print("Atomic coordinates:", coord)
    print("Atomic charges:", nch)
    nu = n - noc
    return nirrep, n, ej, noc, nu, n2el, natms, coord, nch


def prefock(n, n2el):  # procedure to solve for S^-1/2 matrix
    oneelint = open("oneelint", "r")
    twoelint = open("twoelnew", "r")
    dipole = open("dipole", "r")
    ovrl_raw = np.fromfile(oneelint, dtype="float64", count=n * n)
    ekin_raw = np.fromfile(oneelint, dtype="float64", count=n * n)
    epot_raw = np.fromfile(oneelint, dtype="float64", count=n * n)
    dipole_raw = np.fromfile(dipole, dtype="float64", count=(-1))
    ovrl = np.reshape(ovrl_raw, (n, n))  # S matrix
    ekin = np.reshape(ekin_raw, (n, n))  # Ek matrix
    epot = np.reshape(epot_raw, (n, n))  # Ep matrix

    ialone = np.int32(65535)
    ibitwd = np.int32(16)
    buf = np.array([])
    ibuf = np.array([])
    nint = np.array([0])
    nint2 = np.array([])
    while nint[0] == 0 or nint[-1] == 600:
        if nint[0] == 0:
            nint = np.array([])
        x = np.fromfile(twoelint, dtype="float64", count=600)
        buf = np.append(buf, x)
        x = np.fromfile(twoelint, dtype="int32", count=1200)
        ibuf = np.append(ibuf, x)
        x = np.fromfile(twoelint, dtype="int32", count=1)
        nint = np.append(nint, x)
        x = np.fromfile(twoelint, dtype="int32", count=1)
        nint2 = np.append(nint2, x)

    # Vectorized index decoding
    ibuf_int = ibuf.astype(np.int32)
    m1_all = ibuf_int[0::2][:n2el]
    m2_all = ibuf_int[1::2][:n2el]
    idx_i = (m1_all & ialone) - 1
    idx_j = ((m1_all >> ibitwd) & ialone) - 1
    idx_k = (m2_all & ialone) - 1
    idx_l = ((m2_all >> ibitwd) & ialone) - 1
    buf_vals = buf[:n2el]

    with open("twoeltemp", "wb") as twoelfile:
        for im in range(0, n):
            twoel_raw = np.zeros((n, n, n))
            # Vectorized: select integrals where l == im
            mask_l = (idx_l == im)
            if np.any(mask_l):
                twoel_raw[idx_i[mask_l], idx_j[mask_l], idx_k[mask_l]] = buf_vals[mask_l]
                twoel_raw[idx_j[mask_l], idx_i[mask_l], idx_k[mask_l]] = buf_vals[mask_l]
            # Vectorized: select integrals where k == im
            mask_k = (idx_k == im)
            if np.any(mask_k):
                twoel_raw[idx_i[mask_k], idx_j[mask_k], idx_l[mask_k]] = buf_vals[mask_k]
                twoel_raw[idx_j[mask_k], idx_i[mask_k], idx_l[mask_k]] = buf_vals[mask_k]
            # Vectorized: select integrals where j == im
            mask_j = (idx_j == im)
            if np.any(mask_j):
                twoel_raw[idx_k[mask_j], idx_l[mask_j], idx_i[mask_j]] = buf_vals[mask_j]
                twoel_raw[idx_l[mask_j], idx_k[mask_j], idx_i[mask_j]] = buf_vals[mask_j]
            # Vectorized: select integrals where i == im
            mask_i = (idx_i == im)
            if np.any(mask_i):
                twoel_raw[idx_k[mask_i], idx_l[mask_i], idx_j[mask_i]] = buf_vals[mask_i]
                twoel_raw[idx_l[mask_i], idx_k[mask_i], idx_j[mask_i]] = buf_vals[mask_i]
            twoel_raw.tofile(twoelfile)  # output of two-electron integrals into the file
    twoelfile.close()

    d, u = np.linalg.eig(ovrl)  # D matrix, U matrix
    d12 = np.transpose(u) @ ovrl @ u
    d12 = np.diag(np.diagonal(d12) ** -0.5)  # D^-1/2 matrix
    s12 = u @ d12 @ np.transpose(u)  # S^-1/2 matrix

    return s12, ovrl, ekin, epot, dipole_raw


def fock(ej, n, ekin, epot, noc, s12):  # RHF main procedure
    niter = 100
    thresh = 1.0e-7
    damping = 0.5  # damping parameter to enhance convergence
    escf = ej
    ca = np.zeros((n, n))  # SCF coefficients matrix
    pa = np.zeros((n, n))
    h = ekin + epot  # single-electron operator matrix E_kin + E_pot
    escfold = 0
    iteration = 1
    fap = np.zeros(n)

    print("Commencing RHF calculations...")

    # Cache two-electron integrals in memory to avoid repeated disk I/O
    twoel_cache = []
    twoelfile = open("twoeltemp", "rb")
    twoelfile.seek(0)
    for l in range(0, n):
        twoel_raw = np.fromfile(twoelfile, dtype="float64", count=n * n * n)
        twoel_cache.append(np.reshape(twoel_raw, (n, n, n)))
    twoelfile.close()

    while np.absolute(escf - escfold) > thresh and iteration < niter:
        if iteration != 1:
            escfold = escf

        pa_new = 2.0 * ca[:, :noc] @ ca[:, :noc].T
        pa = damping * pa + (1 - damping) * pa_new

        x = np.zeros((n, n))
        for l in range(0, n):
            twoel = twoel_cache[l]
            x += np.einsum('k,ijk->ij', pa[:, l], twoel)
            x[:, l] += -0.5 * np.einsum('kj,ijk->i', pa, twoel)

        fa = h + x  # Fock matrix

        escf = 0.5 * pa * (h + fa)
        escf = ej + np.sum(escf)

        print("Iteration =", iteration, "      E_SCF =", escf)

        fap_raw = np.transpose(s12) @ fa @ s12
        fap, cfa = np.linalg.eig(fap_raw)  # fap = orbital energies, cfa = eigenvectors F'alpha
        idx = fap.argsort()
        fap = fap[idx]
        cfa = cfa[:, idx]
        ca = s12 @ cfa

        iteration = iteration + 1

    if iteration == niter:
        print("Couldn't converge in", iteration, "iterations!")
        print("Aborting execution of the program!")
        sys.exit("Did NOT converge RHF")
    else:
        print("RHF energy =", escf)
        print("Converged in", iteration - 1, "iterations")
    return ca, escf, fap, pa


def trans4(n, ca):  # four-index transformation, used in post-HF methods
    temp = np.zeros((n, n, n, n))

    # Cache two-electron integrals in memory — avoid re-opening file n times
    twoel_cache = []
    twoelfile = open("twoeltemp", "rb")
    twoelfile.seek(0)
    for u in range(0, n):
        twoel_raw = np.fromfile(twoelfile, dtype="float64", count=n * n * n)
        twoel_cache.append(np.reshape(twoel_raw, (n, n, n)))
    twoelfile.close()

    with open("vmoltemp", "wb") as vmolfile:
        for l in range(0, n):
            for u in range(0, n):
                temp[:, :, :, l] = np.add(temp[:, :, :, l], ca[u, l] * twoel_cache[u], out=temp[:, :, :, l],
                                          casting="unsafe")
                # Step 2: contract over t -> temp2[:,:,k,l] = sum_t ca[t,k] * temp[:,: ,t,l]
            temp2_l = np.einsum('ijt,tk->ijk', temp[:, :, :, l], ca, optimize=True)
            # Step 3: contract over s -> temp3[:,j,k,l] = sum_s ca[s,j] * temp2[:,s,k,l]
            temp3_l = np.einsum('isk,sj->ijk', temp2_l, ca, optimize=True)
            # Step 4: contract over r -> vmol[i,j,k,l] = sum_r ca[r,i] * temp3[r,j,k,l]
            vmol_l = np.einsum('rjk,ri->ijk', temp3_l, ca, optimize=True)
            vmol_l.tofile(vmolfile)  # exporting molecular orbitals into the file

    vmolfile.close()

    return None


def sortab(n, noc, nu):  # sorting

    vmolfile = open("vmoltemp", "rb")

    vmolfile.seek(n * n * n * 8 * noc)
    vp = np.zeros((nu, nu, nu, nu))
    for d in range(0, nu):
        vmol_raw = np.fromfile(vmolfile, dtype="float64", count=n * n * n)
        vmol = np.reshape(vmol_raw, (n, n, n))
        vmol = np.transpose(vmol, (0, 2, 1))
        vp[:, :, :, d] = vmol[noc:, noc:, noc:]

    vmolfile.seek(0)
    vh = np.zeros((noc, noc, noc, noc))
    vl = np.zeros((noc, nu, nu, noc))
    vr = np.zeros((noc, nu, nu, noc))
    for l in range(0, noc):
        vmol_raw = np.fromfile(vmolfile, dtype="float64", count=n * n * n)
        vmol = np.reshape(vmol_raw, (n, n, n))
        vmol = np.transpose(vmol, (0, 2, 1))
        vh[:, :, :, l] = vmol[:noc, :noc, :noc]
        vmol = np.transpose(vmol, (1, 0, 2))
        vl[:, :, :, l] = vmol[:noc, noc:, noc:]
        vmol = np.transpose(vmol, (1, 2, 0))
        vr[:, :, :, l] = vmol[:noc, noc:, noc:]

    vmolfile.close()

    return vh, vp, vl, vr


def mp2(noc, nu, fap, vr, escf):  # MP2
    # Vectorized denominator tensor via broadcasting
    fi = fap[:noc]
    fa = fap[noc:noc + nu]
    den = (fi[:, None, None, None] + fi[None, None, None, :]
           - fa[None, :, None, None] - fa[None, None, :, None])

    mp2 = vr * (2 * vr - np.transpose(vr, (0, 2, 1, 3))) / den
    mp2 = np.sum(mp2)
    c2 = vr / den

    print("MP2 correction is equal to:", mp2)
    print("Total MP2 energy:", escf + mp2)
    return c2, den, mp2


def mp3(c2, vr, vl, vh, vp, escf, mp2):  # MP3
    temp1 = np.einsum('mebj,meai->iabj', c2, vr, optimize=True) * 2 + np.einsum('meai,mebj->iabj', c2, vr,
                                                                                optimize=True) * 2
    temp2 = np.einsum('mbej,meai->iabj', c2, vr, optimize=True) * -1 - np.einsum('maei,mebj->iabj', c2, vr,
                                                                                 optimize=True)
    temp3 = np.einsum('mebj,iaem->iabj', c2, vl, optimize=True) * -1 - np.einsum('meai,jbem->iabj', c2, vl,
                                                                                 optimize=True)
    temp4 = np.einsum('maej,ibem->iabj', c2, vl, optimize=True) * -1 - np.einsum('mbei,jaem->iabj', c2, vl,
                                                                                 optimize=True)
    temp5 = np.einsum('mabn,mnij->iabj', c2, vh, optimize=True)
    temp6 = np.einsum('iefj,efab->iabj', c2, vp, optimize=True)
    c2n = temp1 + temp2 + temp3 + temp4 + temp5 + temp6
    mp3 = c2n * (2 * c2 - np.transpose(c2, (0, 2, 1, 3)))
    mp3 = np.sum(mp3)

    print("MP3 correction is equal to:", mp3)
    print("Total MP3 energy:", escf + mp2 + mp3)
    return mp3


def lccd(c2, mp2, vr, vl, vh, vp, escf, den):
    niter = 100
    thresh = 1.0e-7
    ecorr_lccd = mp2
    c2o = c2
    ecorrold = 0
    iteration = 1
    while np.absolute(ecorr_lccd - ecorrold) > thresh and iteration < niter:
        if iteration != 1:
            ecorrold = ecorr_lccd

        temp1 = np.einsum('mebj,meai->iabj', c2o, vr, optimize=True) * 2
        temp2 = np.einsum('mbej,meai->iabj', c2o, vr, optimize=True) * -1
        temp3 = np.einsum('mebj,iaem->iabj', c2o, vl, optimize=True) * -1
        temp4 = np.einsum('maej,ibem->iabj', c2o, vl, optimize=True) * -1
        temp5 = np.einsum('mabn,mnij->iabj', c2o, vh, optimize=True) * 0.5
        temp6 = np.einsum('iefj,efab->iabj', c2o, vp, optimize=True) * 0.5
        c2n = temp1 + temp2 + temp3 + temp4 + temp5 + temp6
        c2n = (c2n + np.transpose(c2n, (3, 2, 1, 0)) + vr) / den
        ecorr_lccd = c2n * (2 * vr - np.transpose(vr, (0, 2, 1, 3)))
        ecorr_lccd = np.sum(ecorr_lccd)
        c2o = c2n

        print("Iteration =", iteration, "      correlation energy =", ecorr_lccd)

        iteration = iteration + 1

    lccd = escf + ecorr_lccd
    if iteration == niter:
        print("Couldn't converge in", iteration, "iterations!")
        print("Aborting execution of the program!")
        sys.exit("Did NOT converge LCCD")
    else:
        print("LCCD energy =", lccd)
        print("Converged in", iteration - 1, "iterations")
    return lccd, ecorr_lccd


def cid(c2, mp2, vr, vl, vh, vp, escf, den):
    niter = 100
    thresh = 1.0e-7
    ecorr_cid = mp2
    c2o = c2
    ecorrold = 0
    iteration = 1
    while np.absolute(ecorr_cid - ecorrold) > thresh and iteration < niter:
        if iteration != 1:
            ecorrold = ecorr_cid

        temp1 = np.einsum('mebj,meai->iabj', c2o, vr, optimize=True) * 2
        temp2 = np.einsum('mbej,meai->iabj', c2o, vr, optimize=True) * -1
        temp3 = np.einsum('mebj,iaem->iabj', c2o, vl, optimize=True) * -1
        temp4 = np.einsum('maej,ibem->iabj', c2o, vl, optimize=True) * -1
        temp5 = np.einsum('mabn,mnij->iabj', c2o, vh, optimize=True) * 0.5
        temp6 = np.einsum('iefj,efab->iabj', c2o, vp, optimize=True) * 0.5
        c2n = temp1 + temp2 + temp3 + temp4 + temp5 + temp6
        c2n = (c2n + np.transpose(c2n, (3, 2, 1, 0)) + vr) / (den + ecorr_cid)
        ecorr_cid = c2n * (2 * vr - np.transpose(vr, (0, 2, 1, 3)))
        ecorr_cid = np.sum(ecorr_cid)
        c2o = c2n

        print("Iteration =", iteration, "      correlation energy =", ecorr_cid)

        iteration = iteration + 1

    cid = escf + ecorr_cid
    if iteration == niter:
        print("Couldn't converge in", iteration, "iterations!")
        print("Aborting execution of the program!")
        sys.exit("Did NOT converge CID")
    else:
        print("CID energy =", cid)
        print("Converged in", iteration - 1, "iteration")
    return cid, ecorr_cid


def dipole(natms, coord, nch, n, dipole_raw, ca, noc, pa):  # dipole moment
    print("Calculating the RHF dipole moment...")

    nch_arr = np.array(nch)
    dipnu = coord @ nch_arr  # nuclei dipole moment (vectorized)

    t = 0
    dipel = np.zeros(3)
    dipel2 = np.zeros(3)
    for p in range(0, 3):
        datm = np.zeros((n, n))
        # Vectorized triangular fill for datm
        r_idx, s_idx = np.tril_indices(n)
        size = len(r_idx)
        datm[r_idx, s_idx] = dipole_raw[t:t + size]
        datm[s_idx, r_idx] = dipole_raw[t:t + size]
        t += size

        dmol = ca.T @ datm @ ca  # AO->MO transformation (vectorized)

        dipel[p] = 2.0 * np.trace(dmol[:noc, :noc])  # electron dipole moment calculated with method 1 (vectorized)

        dipel2[p] = np.sum(pa * datm)  # electron dipole moment calculated with method 2 (vectorized)

    dipole_moment = np.zeros(3)
    for i in range(0, 3):
        dipole_moment[i] = dipnu[i] + dipel[i]
        if np.absolute(dipole_moment[i]) < 10e-12:
            dipole_moment[i] = 0
    print("Dipole moment, X component:", dipole_moment[0])
    print("Dipole moment, Y component:", dipole_moment[1])
    print("Dipole moment, Z component:", dipole_moment[2])
    total_dipole_moment = np.sqrt(dipole_moment[0] ** 2 + dipole_moment[1] ** 2 + dipole_moment[2] ** 2)
    print("Total dipole moment:", total_dipole_moment)
    return dipole_moment, total_dipole_moment, dipnu, dipel


def quadrupole(natms, nch, coord, n, dipole_raw, pa):  # quadrupole moment
    print("Calculating the RHF quadruple moment...")

    tempquad = np.zeros(6)
    quadnuc = np.zeros(6)
    a = 0
    for ic in range(0, 3):
        for jc in range(0, ic + 1):
            for i in range(0, natms):
                rxx = nch[i] * coord[ic, i] * coord[jc, i]
                tempquad[a] += rxx
            a += 1
    xx = (tempquad[0] + tempquad[2] + tempquad[5]) / 2
    for i in range(0, 6):
        quadnuc[i] = 1.5 * tempquad[i]
    quadnuc[0] = quadnuc[0] - xx
    quadnuc[2] = quadnuc[2] - xx
    quadnuc[5] = quadnuc[5] - xx
    for i in range(0, 6):
        tempquad[i] = quadnuc[i]
    quadnuc[0] = tempquad[0]
    quadnuc[1] = tempquad[2]
    quadnuc[2] = tempquad[5]
    quadnuc[3] = tempquad[1]
    quadnuc[4] = tempquad[3]
    quadnuc[5] = tempquad[4]  # Nuclei quadruple moment

    t = n * (n + 1) // 2 * 3  # vectorized offset calculation

    quadel = np.zeros(6)
    for p in range(0, 6):
        datm = np.zeros((n, n))
        # Vectorized triangular fill for datm
        r_idx, s_idx = np.tril_indices(n)
        size = len(r_idx)
        datm[r_idx, s_idx] = dipole_raw[t:t + size]
        datm[s_idx, r_idx] = dipole_raw[t:t + size]
        t += size
        quadel[p] = np.sum(pa * datm)  # vectorized quadrupole electron contribution

    quadrupole_moment = np.zeros(6)
    for i in range(0, 6):
        quadrupole_moment[i] = quadnuc[i] + quadel[i]
        if np.absolute(quadrupole_moment[i]) < 10e-12:
            quadrupole_moment[i] = 0
    print("Quadruple moment, XX component:", quadrupole_moment[0])
    print("Quadruple moment, YY component:", quadrupole_moment[1])
    print("Quadruple moment, ZZ component:", quadrupole_moment[2])
    print("Quadruple moment, XY component:", quadrupole_moment[3])
    print("Quadruple moment, XZ component:", quadrupole_moment[4])
    print("Quadruple moment, YZ component:", quadrupole_moment[5])
    return quadrupole_moment


def finite_field(n, dipole_raw, ekin, epot, ej, dipnu, noc, s12, nu, vr, vl, vh, vp, escf, mp2,
                 mp3):  # finite-field approach
    field = 0.0002
    dipole_moment_ff_hf = [0, 0, 0]
    alpha_ff_hf = [0, 0, 0]
    dipole_moment_ff_mp2 = [0, 0, 0]
    alpha_ff_mp2 = [0, 0, 0]
    dipole_moment_ff_mp3 = [0, 0, 0]
    alpha_ff_mp3 = [0, 0, 0]
    hfield = np.zeros((n, n))
    escf_x = []
    escf_y = []
    escf_z = []
    mp2_x = []
    mp2_y = []
    mp2_z = []
    mp3_x = []
    mp3_y = []
    mp3_z = []

    datm = np.zeros((3, n, n))
    t = 0
    for p in range(0, 3):
        for r in range(0, n):
            for s in range(0, r + 1):
                datm[p, r, s] = dipole_raw[t]
                datm[p, s, r] = dipole_raw[t]
                t = t + 1

    # Cache two-electron integrals ONCE for all finite-field SCF iterations
    twoel_cache_ff = []
    twoelfile = open("twoeltemp", "rb")
    twoelfile.seek(0)
    for l in range(0, n):
        twoel_raw = np.fromfile(twoelfile, dtype="float64", count=n * n * n)
        twoel_cache_ff.append(np.reshape(twoel_raw, (n, n, n)))
    twoelfile.close()

    for a in range(0, 3):
        if a == 0:
            print("Calculating finite field energy for the X direction...")
        if a == 1:
            print("Calculating finite field energy for the Y direction...")
        if a == 2:
            print("Calculating finite field energy for the Z direction...")
        for b in range(0, 4):
            if b == 0:
                efield = field
            elif b == 1:
                efield = field * -1
            elif b == 2:
                efield = field * 2
            else:
                efield = field * -2

            niter = 200
            thresh = 1.0e-7
            damping = 0.5
            ca_new = np.zeros((n, n))  # SCF coefficients matrix
            pa_new = np.zeros((n, n))
            h = ekin + epot  # single-electron operator matrix E_kin + E_pot
            escfold = 0
            iteration = 1
            hfield = h + efield * datm[a]  # vectorized field Hamiltonian

            ejfield = ej + efield * dipnu[a]
            escf_new = ejfield
            fap_new = np.zeros(n)

            while np.absolute(escf_new - escfold) > thresh and iteration < niter:
                if iteration != 1:
                    escfold = escf_new

                pa_calc = 2.0 * ca_new[:, :noc] @ ca_new[:, :noc].T
                pa_new = damping * pa_new + (1 - damping) * pa_calc

                x = np.zeros((n, n))
                for l in range(0, n):
                    twoel = twoel_cache_ff[l]
                    x += np.einsum('k,ijk->ij', pa_new[:, l], twoel)
                    x[:, l] += -0.5 * np.einsum('kj,ijk->i', pa_new, twoel)

                fa_new = hfield + x  # Fock matrix

                escf_new = 0.5 * pa_new * (hfield + fa_new)
                escf_new = ejfield + np.sum(escf_new)

                fap_raw = np.transpose(s12) @ fa_new @ s12
                fap_new, cfa_new = np.linalg.eig(fap_raw)  # fap = orbital energies, cfa = eigenvectors F'alpha
                idx = fap_new.argsort()
                fap_new = fap_new[idx]
                cfa_new = cfa_new[:, idx]
                ca_new = s12 @ cfa_new

                iteration = iteration + 1

            if a == 0:
                escf_x.append(escf_new)
            if a == 1:
                escf_y.append(escf_new)
            if a == 2:
                escf_z.append(escf_new)

            # Vectorized denominator tensor via broadcasting
            fi_ff = fap_new[:noc]
            fa_ff = fap_new[noc:noc + nu]
            den = (fi_ff[:, None, None, None] + fi_ff[None, None, None, :]
                   - fa_ff[None, :, None, None] - fa_ff[None, None, :, None])

            mp2_new = vr * (2 * vr - np.transpose(vr, (0, 2, 1, 3))) / den
            mp2_new = np.sum(mp2_new)
            c2_new = vr / den

            if a == 0:
                mp2_x.append(escf_new + mp2_new)
            if a == 1:
                mp2_y.append(escf_new + mp2_new)
            if a == 2:
                mp2_z.append(escf_new + mp2_new)

            temp1 = np.einsum('mebj,meai->iabj', c2_new, vr, optimize=True) * 2 + np.einsum('meai,mebj->iabj', c2_new,
                                                                                            vr, optimize=True) * 2
            temp2 = np.einsum('mbej,meai->iabj', c2_new, vr, optimize=True) * -1 - np.einsum('maei,mebj->iabj', c2_new,
                                                                                             vr, optimize=True)
            temp3 = np.einsum('mebj,iaem->iabj', c2_new, vl, optimize=True) * -1 - np.einsum('meai,jbem->iabj', c2_new,
                                                                                             vl, optimize=True)
            temp4 = np.einsum('maej,ibem->iabj', c2_new, vl, optimize=True) * -1 - np.einsum('mbei,jaem->iabj', c2_new,
                                                                                             vl, optimize=True)
            temp5 = np.einsum('mabn,mnij->iabj', c2_new, vh, optimize=True)
            temp6 = np.einsum('iefj,efab->iabj', c2_new, vp, optimize=True)
            c2n_new = temp1 + temp2 + temp3 + temp4 + temp5 + temp6
            mp3_new = c2n_new * (2 * c2_new - np.transpose(c2_new, (0, 2, 1, 3)))
            mp3_new = np.sum(mp3_new)

            if a == 0:
                mp3_x.append(escf_new + mp2_new + mp3_new)
            if a == 1:
                mp3_y.append(escf_new + mp2_new + mp3_new)
            if a == 2:
                mp3_z.append(escf_new + mp2_new + mp3_new)

    # For the HF method
    dipole_moment_ff_hf[0] = -1 * (8 * (escf_x[0] - escf_x[1]) - (escf_x[2] - escf_x[3])) / (12 * field)
    dipole_moment_ff_hf[1] = -1 * (8 * (escf_y[0] - escf_y[1]) - (escf_y[2] - escf_y[3])) / (12 * field)
    dipole_moment_ff_hf[2] = -1 * (8 * (escf_z[0] - escf_z[1]) - (escf_z[2] - escf_z[3])) / (12 * field)
    total_dipole_moment_ff_hf = np.sqrt(
        dipole_moment_ff_hf[0] ** 2 + dipole_moment_ff_hf[1] ** 2 + dipole_moment_ff_hf[2] ** 2)

    for i in range(0, 3):
        if np.absolute(dipole_moment_ff_hf[i]) < 10e-8:
            dipole_moment_ff_hf[i] = 0

    alpha_ff_hf[0] = (
                -1 * (((16 * (escf_x[0] + escf_x[1])) - (escf_x[2] + escf_x[3]) - (30 * escf)) / (12 * field ** 2)))
    alpha_ff_hf[1] = (
                -1 * (((16 * (escf_y[0] + escf_y[1])) - (escf_y[2] + escf_y[3]) - (30 * escf)) / (12 * field ** 2)))
    alpha_ff_hf[2] = (
                -1 * (((16 * (escf_z[0] + escf_z[1])) - (escf_z[2] + escf_z[3]) - (30 * escf)) / (12 * field ** 2)))
    alpha_ave_ff_hf = (alpha_ff_hf[0] + alpha_ff_hf[1] + alpha_ff_hf[2]) / 3

    print("RHF electric properties obtained with the finite-field approach:")
    print("Dipole moment, X component:", dipole_moment_ff_hf[0])
    print("Dipole moment, Y component:", dipole_moment_ff_hf[1])
    print("Dipole moment, Z component:", dipole_moment_ff_hf[2])
    print("Total dipole moment:", total_dipole_moment_ff_hf)
    print("Polarizability, XX component:", alpha_ff_hf[0])
    print("Polarizability, YY component:", alpha_ff_hf[1])
    print("Polarizability, ZZ component:", alpha_ff_hf[2])
    print("Average polarizability:", alpha_ave_ff_hf)

    # For the MP2 method
    dipole_moment_ff_mp2[0] = -1 * (8 * (mp2_x[0] - mp2_x[1]) - (mp2_x[2] - mp2_x[3])) / (12 * field)
    dipole_moment_ff_mp2[1] = -1 * (8 * (mp2_y[0] - mp2_y[1]) - (mp2_y[2] - mp2_y[3])) / (12 * field)
    dipole_moment_ff_mp2[2] = -1 * (8 * (mp2_z[0] - mp2_z[1]) - (mp2_z[2] - mp2_z[3])) / (12 * field)
    total_dipole_moment_ff_mp2 = np.sqrt(
        dipole_moment_ff_mp2[0] ** 2 + dipole_moment_ff_mp2[1] ** 2 + dipole_moment_ff_mp2[2] ** 2)

    for i in range(0, 3):
        if np.absolute(dipole_moment_ff_mp2[i]) < 10e-8:
            dipole_moment_ff_mp2[i] = 0

    alpha_ff_mp2[0] = (
                -1 * (((16 * (mp2_x[0] + mp2_x[1])) - (mp2_x[2] + mp2_x[3]) - (30 * (escf + mp2))) / (12 * field ** 2)))
    alpha_ff_mp2[1] = (
                -1 * (((16 * (mp2_y[0] + mp2_y[1])) - (mp2_y[2] + mp2_y[3]) - (30 * (escf + mp2))) / (12 * field ** 2)))
    alpha_ff_mp2[2] = (
                -1 * (((16 * (mp2_z[0] + mp2_z[1])) - (mp2_z[2] + mp2_z[3]) - (30 * (escf + mp2))) / (12 * field ** 2)))
    alpha_ave_ff_mp2 = (alpha_ff_mp2[0] + alpha_ff_mp2[1] + alpha_ff_mp2[2]) / 3

    print("MP2 electric properties obtained with the finite-field approach:")
    print("Dipole moment, X component:", dipole_moment_ff_mp2[0])
    print("Dipole moment, Y component:", dipole_moment_ff_mp2[1])
    print("Dipole moment, Z component:", dipole_moment_ff_mp2[2])
    print("Total dipole moment:", total_dipole_moment_ff_mp2)
    print("Polarizability, XX component:", alpha_ff_mp2[0])
    print("Polarizability, YY component:", alpha_ff_mp2[1])
    print("Polarizability, ZZ component:", alpha_ff_mp2[2])
    print("Average polarizability:", alpha_ave_ff_mp2)

    # For the MP3 method
    dipole_moment_ff_mp3[0] = -1 * (8 * (mp3_x[0] - mp3_x[1]) - (mp3_x[2] - mp3_x[3])) / (12 * field)
    dipole_moment_ff_mp3[1] = -1 * (8 * (mp3_y[0] - mp3_y[1]) - (mp3_y[2] - mp3_y[3])) / (12 * field)
    dipole_moment_ff_mp3[2] = -1 * (8 * (mp3_z[0] - mp3_z[1]) - (mp3_z[2] - mp3_z[3])) / (12 * field)
    total_dipole_moment_ff_mp3 = np.sqrt(
        dipole_moment_ff_mp3[0] ** 2 + dipole_moment_ff_mp3[1] ** 2 + dipole_moment_ff_mp3[2] ** 2)

    for i in range(0, 3):
        if np.absolute(dipole_moment_ff_mp3[i]) < 10e-8:
            dipole_moment_ff_mp3[i] = 0

    alpha_ff_mp3[0] = (-1 * (
                ((16 * (mp3_x[0] + mp3_x[1])) - (mp3_x[2] + mp3_x[3]) - (30 * (escf + mp2 + mp3))) / (12 * field ** 2)))
    alpha_ff_mp3[1] = (-1 * (
                ((16 * (mp3_y[0] + mp3_y[1])) - (mp3_y[2] + mp3_y[3]) - (30 * (escf + mp2 + mp3))) / (12 * field ** 2)))
    alpha_ff_mp3[2] = (-1 * (
                ((16 * (mp3_z[0] + mp3_z[1])) - (mp3_z[2] + mp3_z[3]) - (30 * (escf + mp2 + mp3))) / (12 * field ** 2)))
    alpha_ave_ff_mp3 = (alpha_ff_mp3[0] + alpha_ff_mp3[1] + alpha_ff_mp3[2]) / 3

    print("MP3 electric properties obtained with the finite-field approach:")
    print("Dipole moment, X component:", dipole_moment_ff_mp3[0])
    print("Dipole moment, Y component:", dipole_moment_ff_mp3[1])
    print("Dipole moment, Z component:", dipole_moment_ff_mp3[2])
    print("Total dipole moment:", total_dipole_moment_ff_mp3)
    print("Polarizability, XX component:", alpha_ff_mp3[0])
    print("Polarizability, YY component:", alpha_ff_mp3[1])
    print("Polarizability, ZZ component:", alpha_ff_mp3[2])
    print("Average polarizability:", alpha_ave_ff_mp3)

    return dipole_moment_ff_hf, total_dipole_moment_ff_hf, alpha_ff_hf, alpha_ave_ff_hf, dipole_moment_ff_mp2, total_dipole_moment_ff_mp2, alpha_ff_mp2, alpha_ave_ff_mp2, dipole_moment_ff_mp3, total_dipole_moment_ff_mp3, alpha_ff_mp3, alpha_ave_ff_mp3
