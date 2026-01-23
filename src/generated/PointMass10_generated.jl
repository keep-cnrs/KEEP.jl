
"""
/!\\ dalpha should ALWAYS be set to 0.
It is not used in the calculations, just returned as is. /!\\

Copy-pasted code from the old code, to be replaced in later iterations.
"""
function vitesse_init(alpha, phi2, theta2, R, tau, dalpha, dtau, p)
	(; r) = p
	theta_0, phi_0, delta_theta, delta_phi = p.θ0, p.φ0, p.Δθ, p.Δφ
	theta = theta_0 + delta_theta * sin(2*tau)
	phi = phi_0 + delta_phi * sin(tau)
	dtheta = 2 *dtau* delta_theta * cos(2*tau)
	dphi =  dtau* delta_phi*cos(tau)
	dphi2, dtheta2, dR = SA[
		-cos(theta2)*sin(theta)^2/r*(sin(phi2+alpha)*cos(phi)-sin(phi)*cos(phi2+alpha))/sin(theta2)/(sin(phi2+alpha)^2*cos(theta2)*cos(theta)+sin(phi2+alpha)*sin(phi)*sin(theta)*sin(theta2)+cos(theta2)*cos(phi2+alpha)^2*cos(theta)+cos(phi2+alpha)*cos(phi)*sin(theta)*sin(theta2))*R*dtheta-0.1e1/r*(cos(theta2)*cos(phi2+alpha)*cos(theta)+cos(phi)*sin(theta)*sin(theta2))/sin(theta2)/(sin(phi2+alpha)^2*cos(theta2)*cos(theta)+sin(phi2+alpha)*sin(phi)*sin(theta)*sin(theta2)+cos(theta2)*cos(phi2+alpha)^2*cos(theta)+cos(phi2+alpha)*cos(phi)*sin(theta)*sin(theta2))*(-R*dtheta*sin(phi)*cos(theta)-R*dphi*cos(phi)*sin(theta))+0.1e1/r*(sin(phi2+alpha)*cos(theta2)*cos(theta)+sin(phi)*sin(theta)*sin(theta2))/sin(theta2)/(sin(phi2+alpha)^2*cos(theta2)*cos(theta)+sin(phi2+alpha)*sin(phi)*sin(theta)*sin(theta2)+cos(theta2)*cos(phi2+alpha)^2*cos(theta)+cos(phi2+alpha)*cos(phi)*sin(theta)*sin(theta2))*(R*dphi*sin(phi)*sin(theta)-R*dtheta*cos(phi)*cos(theta))
		sin(theta)^2/r*(sin(phi2+alpha)*sin(phi)+cos(phi2+alpha)*cos(phi))/(sin(phi2+alpha)^2*cos(theta2)*cos(theta)+sin(phi2+alpha)*sin(phi)*sin(theta)*sin(theta2)+cos(theta2)*cos(phi2+alpha)^2*cos(theta)+cos(phi2+alpha)*cos(phi)*sin(theta)*sin(theta2))*R*dtheta-cos(theta)*sin(phi2+alpha)/r/(sin(phi2+alpha)^2*cos(theta2)*cos(theta)+sin(phi2+alpha)*sin(phi)*sin(theta)*sin(theta2)+cos(theta2)*cos(phi2+alpha)^2*cos(theta)+cos(phi2+alpha)*cos(phi)*sin(theta)*sin(theta2))*(-R*dtheta*sin(phi)*cos(theta)-R*dphi*cos(phi)*sin(theta))-cos(theta)*cos(phi2+alpha)/r/(sin(phi2+alpha)^2*cos(theta2)*cos(theta)+sin(phi2+alpha)*sin(phi)*sin(theta)*sin(theta2)+cos(theta2)*cos(phi2+alpha)^2*cos(theta)+cos(phi2+alpha)*cos(phi)*sin(theta)*sin(theta2))*(R*dphi*sin(phi)*sin(theta)-R*dtheta*cos(phi)*cos(theta))
		cos(theta2)*(sin(phi2+alpha)^2+cos(phi2+alpha)^2)/(sin(phi2+alpha)^2*cos(theta2)*cos(theta)+sin(phi2+alpha)*sin(phi)*sin(theta)*sin(theta2)+cos(theta2)*cos(phi2+alpha)^2*cos(theta)+cos(phi2+alpha)*cos(phi)*sin(theta)*sin(theta2))*R*dtheta*sin(theta)+sin(theta2)*sin(phi2+alpha)/(sin(phi2+alpha)^2*cos(theta2)*cos(theta)+sin(phi2+alpha)*sin(phi)*sin(theta)*sin(theta2)+cos(theta2)*cos(phi2+alpha)^2*cos(theta)+cos(phi2+alpha)*cos(phi)*sin(theta)*sin(theta2))*(-R*dtheta*sin(phi)*cos(theta)-R*dphi*cos(phi)*sin(theta))+sin(theta2)*cos(phi2+alpha)/(sin(phi2+alpha)^2*cos(theta2)*cos(theta)+sin(phi2+alpha)*sin(phi)*sin(theta)*sin(theta2)+cos(theta2)*cos(phi2+alpha)^2*cos(theta)+cos(phi2+alpha)*cos(phi)*sin(theta)*sin(theta2))*(R*dphi*sin(phi)*sin(theta)-R*dtheta*cos(phi)*cos(theta))
]

	
	return dalpha, dphi2, dtheta2, dR, dtau
end
