function displacementLeft(t)
	return math.sin(2.0 * math.pi * t)
end

function velocityLeft(t)
	return 2.0 * math.pi * math.cos(2.0 * math.pi * t)
end

function displacementRight(t)
	return math.sin(2.0 * math.pi * (t - 1.0))
end

function velocityRight(t)
	return 2.0 * math.pi * math.cos(2.0 * math.pi * (t - 1.0))
end

function initialVel(x)
	return 2.0 * math.pi * math.cos(2.0 * math.pi * x)
end

function initialDisp(x)
	return -math.sin(2.0 * math.pi * x)
end

function bodyForce(x, t)
	local E = 2.0
	local r = 1.0
	local b0 = (E - r) / r
	local two_pi = 2.0 * math.pi
	local four_pi_sqr = 4.0 * math.pi * math.pi
	return b0 * four_pi_sqr * math.sin(two_pi * (t - x))
end

function Ricker_wavelet(x)
	local fp = 1.0
	local amp = 1.0
	local xs = 12.5
	return amp * (1.0 - 2.0 * (math.pi * fp * (x - xs)) ^ 2) * math.exp(-(math.pi * fp * (x - xs)) ^ 2)
end

function exactDisp(x, t)
	local two_pi = 2.0 * math.pi
	return math.sin(two_pi * (t - x))
end

function exactVel(x, t)
	local two_pi = 2.0 * math.pi
	return two_pi * math.cos(two_pi * (t - x))
end
