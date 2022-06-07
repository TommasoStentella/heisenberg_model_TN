function gauge(psi, N, bonds, bond)
    up = prod(psi)
    n = up*up
    n = n[1]
    
    up = up/sqrt(n)
    block = ITensor()
    for i=1:N-1
        ind = unioninds(up)[i]
        if i == 1
            block = up
            U,S,V = svd(block, ind; maxdim=bond)
            psi[i] = U
            bonds[i] = commoninds(U,S)[1]
            block = S*V
        elseif i == N-1
            U,S,V = svd(block, ind, bonds[i-1]; maxdim=bond)
            psi[i] = U
            bonds[i] = commoninds(U,S)[1]
            psi[i+1] = S*V
        else
            U,S,V = svd(block, ind, bonds[i-1]; maxdim=bond)
            psi[i] = U
            bonds[i] = commoninds(U,S)[1]
            block = S*V
        end
    end
    
    return psi, bonds
end


function HXXZ_dmrg(N, Jxx, Jz, sweeps=3, bond=4)
    
    if iseven(N) == false
        N -= 1
    end
    
    o_ind = [Index(2) for i=1:N]
    bonds = [Index(bond) for i=1:N-1]
    
    Sz = ITensor[]
    Sp = ITensor[]
    Sm = ITensor[]
    Id = ITensor[]
    psi = ITensor[]
    psi_bar = ITensor[]
    for i=1:N
        push!(Sz, ITensor(o_ind[i], o_ind[i]'))
        Sz[i][o_ind[i]=>1, o_ind[i]'=>1] = 1
        Sz[i][o_ind[i]=>1, o_ind[i]'=>2] = 0
        Sz[i][o_ind[i]=>2, o_ind[i]'=>1] = 0
        Sz[i][o_ind[i]=>2, o_ind[i]'=>2] = -1
        
        push!(Sp, ITensor(o_ind[i], o_ind[i]'))
        Sp[i][o_ind[i]=>1, o_ind[i]'=>1] = 0
        Sp[i][o_ind[i]=>1, o_ind[i]'=>2] = 1
        Sp[i][o_ind[i]=>2, o_ind[i]'=>1] = 0
        Sp[i][o_ind[i]=>2, o_ind[i]'=>2] = 0
        
        push!(Sm, ITensor(o_ind[i], o_ind[i]'))
        Sm[i][o_ind[i]=>1, o_ind[i]'=>1] = 0
        Sm[i][o_ind[i]=>1, o_ind[i]'=>2] = 1
        Sm[i][o_ind[i]=>2, o_ind[i]'=>1] = 0
        Sm[i][o_ind[i]=>2, o_ind[i]'=>2] = 0
        
        push!(Id, ITensor(o_ind[i], o_ind[i]'))
        Id[i][o_ind[i]=>1, o_ind[i]'=>1] = 1
        Id[i][o_ind[i]=>1, o_ind[i]'=>2] = 0
        Id[i][o_ind[i]=>2, o_ind[i]'=>1] = 0
        Id[i][o_ind[i]=>2, o_ind[i]'=>2] = 1
        
        if i == 1
            push!(psi, randomITensor(o_ind[i], bonds[i]))
        elseif i == N
            push!(psi, randomITensor(o_ind[i], bonds[i-1]))
        else
            push!(psi, randomITensor(o_ind[i], bonds[i-1], bonds[i]))
        end
    end
    
    psi, bonds = gauge(psi, N, bonds, bond)
    for i=1:N
        push!(psi_bar, psi[i]')
    end
    
    
    H = ITensor()
    
    for i=1:N-1
        temp = 0.5 * (Jxx*(Sp[i]*Sm[i+1]+Sm[i]*Sp[i+1]) + Jz*Sz[i]Sz[i+1])
        for j=1:i-1
            temp = temp * Id[j]
        end
        for j=i+2:N
            temp = temp * Id[j]
        end
        H += temp
    end
    
        
    """
    SWEEPING
    """
    energy = 0
    for s=1:sweeps
        for i=2:N
            block = ITensor()
            for j=1:N
                if j == 1
                    block = H
                end
                block = block * psi_bar[j]
                if (j<(i-1)) | (j>i)
                    block = block * psi[j]
                end
            end

            if i == 2
                U,S,V = svd(block, o_ind[1]; maxdim=bond)
                psi[1] = U
                bonds[1] = commoninds(U,S)[1]
                psi_bar[1] = prime(psi[1])
                psi[2] = S*V
                psi_bar[2] = prime(psi[2])
            else
                U,S,V = svd(block, o_ind[i-1], bonds[i-2]; maxdim=bond)
                psi[i-1] = U
                bonds[i-1] = commoninds(U,S)[1]
                psi_bar[i-1] = prime(psi[i-1])
                psi[i] = S*V
                psi_bar[i] = prime(psi[i])
            end
        end
        psi, bonds = gauge(psi, N, bonds, bond)
        for i=1:N
            psi_bar[i] = psi[i]'
        end
        for i=N:-1:2
            block = ITensor()
            for j=1:N
                if j == 1
                    block = H
                end
                block = block * psi_bar[j]
                if (j<(i-1)) | (j>i)
                    block = block * psi[j]
                end
            end
        
            if i == 2
                U,S,V = svd(block, o_ind[1]; maxdim=bond)
                psi[1] = U
                bonds[1] = commoninds(U,S)[1]
                psi_bar[1] = prime(psi[1])
                psi[2] = S*V
                psi_bar[2] = prime(psi[2])
            else
                U,S,V = svd(block, o_ind[i-1], bonds[i-2]; maxdim=bond)
                psi[i-1] = U
                bonds[i-1] = commoninds(U,S)[1]
                psi_bar[i-1] = prime(psi[i-1])
                psi[i] = S*V
                psi_bar[i] = prime(psi[i])
            end
        end
        psi, bonds = gauge(psi, N, bonds, bond)
        for i=1:N
            psi_bar[i] = psi[i]'
        end
        
        energy = H
        for i=1:N
            energy = energy * psi_bar[i] * psi[i]
        end
        print("energy sweep $s: ",energy[1],"\n")
    end
    
    return energy[1]
end


function HXXZ_dmrg_opt(N, Jxx, Jz, sweeps=3, bond=4)
    
    if iseven(N) == false
        N -= 1
    end
    
    o_ind = [Index(2) for i=1:N]
    bonds = [Index(bond) for i=1:N-1]
    
    Sz = ITensor[]
    Sp = ITensor[]
    Sm = ITensor[]
    Id = ITensor[]
    psi = ITensor[]
    psi_bar = ITensor[]
    for i=1:N
        push!(Sz, ITensor(o_ind[i], o_ind[i]'))
        Sz[i][o_ind[i]=>1, o_ind[i]'=>1] = 1
        Sz[i][o_ind[i]=>1, o_ind[i]'=>2] = 0
        Sz[i][o_ind[i]=>2, o_ind[i]'=>1] = 0
        Sz[i][o_ind[i]=>2, o_ind[i]'=>2] = -1
        
        push!(Sp, ITensor(o_ind[i], o_ind[i]'))
        Sp[i][o_ind[i]=>1, o_ind[i]'=>1] = 0
        Sp[i][o_ind[i]=>1, o_ind[i]'=>2] = 1
        Sp[i][o_ind[i]=>2, o_ind[i]'=>1] = 0
        Sp[i][o_ind[i]=>2, o_ind[i]'=>2] = 0
        
        push!(Sm, ITensor(o_ind[i], o_ind[i]'))
        Sm[i][o_ind[i]=>1, o_ind[i]'=>1] = 0
        Sm[i][o_ind[i]=>1, o_ind[i]'=>2] = 1
        Sm[i][o_ind[i]=>2, o_ind[i]'=>1] = 0
        Sm[i][o_ind[i]=>2, o_ind[i]'=>2] = 0
        
        push!(Id, ITensor(o_ind[i], o_ind[i]'))
        Id[i][o_ind[i]=>1, o_ind[i]'=>1] = 1
        Id[i][o_ind[i]=>1, o_ind[i]'=>2] = 0
        Id[i][o_ind[i]=>2, o_ind[i]'=>1] = 0
        Id[i][o_ind[i]=>2, o_ind[i]'=>2] = 1
        
        if i == 1
            push!(psi, randomITensor(o_ind[i], bonds[i]))
        elseif i == N
            push!(psi, randomITensor(o_ind[i], bonds[i-1]))
        else
            push!(psi, randomITensor(o_ind[i], bonds[i-1], bonds[i]))
        end
    end
    
    psi, bonds = gauge(psi, N, bonds, bond)
    for i=1:N
        push!(psi_bar, psi[i]')
    end
    
        
    """
    SWEEPING
    """
    energy = 0
    for s=1:sweeps
        for i=2:N
            block = ITensor()
            prods = [ITensor() for k=1:N-1]
            for j=1:N
                for k=1:N-1
                    if (j == k) & (j == 1)
                        prods[k] = 0.5 * (Jxx*(Sp[k]*Sm[k+1]+Sm[k]*Sp[k+1]) + Jz*Sz[k]Sz[k+1])
                    elseif (j == k) & (j != 1)
                        prods[k] = prods[k] * 0.5 * (Jxx*(Sp[k]*Sm[k+1]+Sm[k]*Sp[k+1]) + Jz*Sz[k]Sz[k+1])
                    elseif (j < k) & (j == 1)
                        prods[k] = Id[j]
                    elseif (j < k) | (j > (k+1))
                        prods[k] = prods[k] * Id[j]
                    end
                    prods[k] = prods[k] * psi_bar[j]
                    if (j<(i-1)) | (j>i)
                        prods[k] = prods[k] * psi[j]
                    end
                end
            end
            block = prods[1]
            for k=1:N-1
                block += prods[k]
            end

            if i == 2
                U,S,V = svd(block, o_ind[1]; maxdim=bond)
                psi[1] = U
                bonds[1] = commoninds(U,S)[1]
                psi_bar[1] = prime(psi[1])
                psi[2] = S*V
                psi_bar[2] = prime(psi[2])
            else
                U,S,V = svd(block, o_ind[i-1], bonds[i-2]; maxdim=bond)
                psi[i-1] = U
                bonds[i-1] = commoninds(U,S)[1]
                psi_bar[i-1] = prime(psi[i-1])
                psi[i] = S*V
                psi_bar[i] = prime(psi[i])
            end
        end
        psi, bonds = gauge(psi, N, bonds, bond)
        for i=1:N
            psi_bar[i] = psi[i]'
        end
        for i=N:-1:2
            block = ITensor()
            prods = [ITensor() for k=1:N-1]
            for j=1:N
                for k=1:N-1
                    if (j == k) & (j == 1)
                        prods[k] = 0.5 * (Jxx*(Sp[k]*Sm[k+1]+Sm[k]*Sp[k+1]) + Jz*Sz[k]Sz[k+1])
                    elseif (j == k) & (j != 1)
                        prods[k] = prods[k] * 0.5 * (Jxx*(Sp[k]*Sm[k+1]+Sm[k]*Sp[k+1]) + Jz*Sz[k]Sz[k+1])
                    elseif (j < k) & (j == 1)
                        prods[k] = Id[j]
                    elseif (j < k) | (j > (k+1))
                        prods[k] = prods[k] * Id[j]
                    end
                    prods[k] = prods[k] * psi_bar[j]
                    if (j<(i-1)) | (j>i)
                        prods[k] = prods[k] * psi[j]
                    end
                end
            end
            block = prods[1]
            for k=1:N-1
                block += prods[k]
            end
        
            if i == 2
                U,S,V = svd(block, o_ind[1]; maxdim=bond)
                psi[1] = U
                bonds[1] = commoninds(U,S)[1]
                psi_bar[1] = prime(psi[1])
                psi[2] = S*V
                psi_bar[2] = prime(psi[2])
            else
                U,S,V = svd(block, o_ind[i-1], bonds[i-2]; maxdim=bond)
                psi[i-1] = U
                bonds[i-1] = commoninds(U,S)[1]
                psi_bar[i-1] = prime(psi[i-1])
                psi[i] = S*V
                psi_bar[i] = prime(psi[i])
            end
        end
        psi, bonds = gauge(psi, N, bonds, bond)
        for i=1:N
            psi_bar[i] = psi[i]'
        end
        
        block = ITensor()
        prods = [ITensor() for k=1:N-1]
        for j=1:N
            for k=1:N-1
                if (j == k) & (j == 1)
                    prods[k] = 0.5 * (Jxx*(Sp[k]*Sm[k+1]+Sm[k]*Sp[k+1]) + Jz*Sz[k]Sz[k+1])
                elseif (j == k) & (j != 1)
                    prods[k] = prods[k] * 0.5 * (Jxx*(Sp[k]*Sm[k+1]+Sm[k]*Sp[k+1]) + Jz*Sz[k]Sz[k+1])
                elseif (j < k) & (j == 1)
                    prods[k] = Id[j]
                elseif (j < k) | (j > (k+1))
                    prods[k] = prods[k] * Id[j]
                end
                prods[k] = prods[k] * psi_bar[j] * psi[j]
            end
        end
        block = prods[1]
        for k=1:N-1
            block += prods[k]
        end
        energy = block[1]
        print("energy sweep $s: ",energy,"\n")
    end
    
    return energy
end