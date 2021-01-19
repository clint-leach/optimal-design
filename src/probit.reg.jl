function probit_reg(y, X, Xfull, nmcmc)

	N = size(Xfull)[1]
	n = length(y)

	beta_save = fill(0.0, 2, nmcmc)
	p_save = fill(0.0, N, nmcmc)
	y_save = BitArray(undef, N, nmcmc)

	# Priors
	Σ_β = ScalMat(2, 2.25)
	Σ_β_inv = inv(Σ_β)
	μ_β = fill(0.0, 2)

	# Initial values
	β = μ_β

	Σtemp = PDMat(Symmetric(inv(Σ_β_inv + X' * X)))

	# Gibbs loop
	for k in 1:nmcmc

		# Sample z
		z = fill(0.0, length(y))
		μ = X * β

		z[y .== 1] = rand.(Truncated.(Normal.(μ[y .== 1], 1.0), 0.0, Inf))
		z[y .== 0] = rand.(Truncated.(Normal.(μ[y .== 0], 1.0), -Inf, 0.0))

		# Sample beta
		βtemp = Σtemp * (Σ_β_inv * μ_β + X' * z)

		β = rand(MvNormal(βtemp, Σtemp))

		beta_save[:, k] = β
		p_save[:, k] = cdf.(Normal(), Xfull * β)
		y_save[:, k] = rand.(Bernoulli.(p_save[:, k]))

	end

	return (β = beta_save, p = p_save, y = y_save)
end


function probit_pprb(ynew, loc, p_prop, β_prop, nmcmc)

	nprop = size(p_prop)[2]

	β_save = fill(0.0, 2, nmcmc)
	y_save = BitArray(undef, N, nmcmc)
	ynew_save = fill(0, nmcmc)

	init = sample(1:nprop, 1)
	p = p_prop[:, init]

	β = β_prop[:, init]
	β_save[:, 1] = β

	@inbounds for k in 1:nmcmc

		idx = sample(1:nprop, 1)
		pstar = p_prop[:, idx]

		mh1 = logpdf(Bernoulli(pstar[loc]), ynew)
		mh2 = logpdf(Bernoulli(p[loc]), ynew)

		r = exp(mh1 - mh2)
		if r > rand(Uniform())
			p = pstar
			β = β_prop[:, idx]
		end

		β_save[:, k] = β
		# p_save[:, k] = p
		y_save[:, k] = rand.(Bernoulli.(p))
		ynew_save[k] = ynew

	end

	return (y = y_save, ynew = ynew_save, β = β_save)

end
