
function designs_pprb(ypred, D, p_prop, β_prop)

	M = size(ypred)[2]
	L = length(D)
	N = size(p_prop)[1]

	# Scores
	vars = fill(0.0, N, L, M)
	means = fill(0.0, N, L, M)

	@progress for j in 1:M
		for l in 1:L

			ynew = ypred[D[l], j]

			updated = probit_pprb(ynew, D[l], p_prop, β_prop, 100000)

			vars[:, l, j] = var(updated.y, dims = 2)
			means[:, l, j] = mean(updated.y, dims = 2)

		end
	end

	return (vars = vars, means = means, score = score)

end
