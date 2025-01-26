
	\subsection{Logarithmic}
	
	Results for Singularity kernel is shown in the figures \ref{fig:Logarithmic_Chol} to \ref{fig:Logarithmic_VCE}. The value of parameter $h$ for this kernel was considered 0.457. Mean and norm ($L_2$) of these three methods are shown in table \ref{tab:Logarithmic_Error}.
	
	\begin{table}[h!]
		\centering
		\caption{Error for interpolation using Logarithmic kernel (unit of values are $\frac{m^2}{s^2}$).}
		\vspace{0.3cm}
		\renewcommand{\arraystretch}{1.4}
		\begin{tabular}{c|c|c|c}
			\textbf{Method} & Cholesky & TSVD & Tikhonov (VCE) \\
			\hline 
			\textbf{Mean Error} & 0.1503 & 0.1746 & 0.4764 \\
			\hline 
			\textbf{Norm of Errors} & 2005.1889 & 2013.2417 & 5436.0035 \\
		\end{tabular}
		\label{tab:Logarithmic_Error}
	\end{table}
	
	\begin{figure}[h!]
		\centering
		\includegraphics[width=16cm]{../Outputs/Logarithmic_Cholesky.pdf}
		\caption{Results of Logarithmic kernel with Cholesky decomposition.}
		\label{fig:Logarithmic_Chol}
	\end{figure}
	
	\begin{figure}[h!]
		\centering
		\includegraphics[width=16cm]{../Outputs/Logarithmic_TSVD.pdf}
		\caption{Results of Logarithmic kernel with TSVD method.}
		\label{fig:Logarithmic_TSVD}
	\end{figure}
	
	\begin{figure}[h!]
		\centering
		\includegraphics[width=16cm]{../Outputs/Logarithmic_VCE.pdf}
		\caption{Results of Logarithmic kernel with Tikhonov (VCE) method.}
		\label{fig:Logarithmic_VCE}
	\end{figure}

	
	\subsection{Abel-Poisson With Noise}
	
	In this case, white-noise with standard deviation of $200 \frac{m^2}{s^2}$ was added to the input values. Results are shown in the figures \ref{fig:AbelPoisson_Chol_Noise} to \ref{fig:AbelPoisson_VCE_Noise}.Also, mean and norm ($L_2$) of the three methods are shown in table \ref{tab:AbelPoisson_Error_Noise}.
	As it can be understood from results, by adding white-noise to the data, Cholesky and TSVD can not achieve satisfactory results. But Tikhonov method using VCE can have much realistic and smoother results that other two methods. 
	
	\begin{table}[h!]
		\centering
		\caption{Error for interpolation using Abel-Poisson kernel and added white-noise (unit of values are $\frac{m^2}{s^2}$).}
		\vspace{0.3cm}
		\renewcommand{\arraystretch}{1.4}
		\begin{tabular}{c|c|c|c}
			\textbf{Method} & Cholesky & TSVD & Tikhonov (VCE) \\
			\hline 
			\textbf{Mean Error} & -3.1558 & -3.1328 & 71.1711 \\
			\hline 
			\textbf{Norm of Errors} & 5548.1711 & 12137.6909 & 8200.7749 \\
		\end{tabular}
		\label{tab:AbelPoisson_Error_Noise}
	\end{table}

	\begin{figure}[h!]
		\centering
		\includegraphics[width=16cm]{../Outputs/AbelPoisson_Cholesky_Noise.pdf}
		\caption{Results of Abel-Poisson kernel with Cholesky decomposition and added white-noise.}
		\label{fig:AbelPoisson_Chol_Noise}
	\end{figure}
	
	\clearpage
	
	\begin{figure}[h!]
		\centering
		\includegraphics[width=16cm]{../Outputs/AbelPoisson_TSVD_Noise.pdf}
		\caption{Results of Abel-Poisson kernel with TSVD method and added white-noise.}
		\label{fig:AbelPoisson_TSVD_Noise}
	\end{figure}
	
	\begin{figure}[h!]
		\centering
		\includegraphics[width=16cm]{../Outputs/AbelPoisson_VCE_Noise.pdf}
		\caption{Results of Abel-Poisson kernel with Tikhonov (VCE) method and added white-noise.}
		\label{fig:AbelPoisson_VCE_Noise}
	\end{figure}
	
	
	\subsection{Singularity With Noise}
	
	The same amount of white-noise was added to the input values. Results are shown in the figures \ref{fig:Singularity_Chol_Noise} to \ref{fig:Singularity_VCE_Noise}.Also, mean and norm ($L_2$) of the three methods are shown in table \ref{tab:Singularity_Error_Noise}.
	
	\begin{table}[h!]
		\centering
		\caption{Error for interpolation using Singularity kernel and added white-noise (unit of values are $\frac{m^2}{s^2}$).}
		\vspace{0.3cm}
		\renewcommand{\arraystretch}{1.4}
		\begin{tabular}{c|c|c|c}
			\textbf{Method} & Cholesky & TSVD & Tikhonov (VCE) \\
			\hline 
			\textbf{Mean Error} & -8.9409 & -8.1949 & -7.2813 \\
			\hline 
			\textbf{Norm of Errors} & 12540.6580 & 5736.3606 & 5770.8221 \\
		\end{tabular}
		\label{tab:Singularity_Error_Noise}
	\end{table}
	
	\begin{figure}[h!]
		\centering
		\includegraphics[width=16cm]{../Outputs/Singularity_Cholesky_Noise.pdf}
		\caption{Results of Singularity kernel with Cholesky decomposition and added white-noise.}
		\label{fig:Singularity_Chol_Noise}
	\end{figure}
	
	\clearpage
	
	\begin{figure}[h!]
		\centering
		\includegraphics[width=16cm]{../Outputs/Singularity_TSVD_Noise.pdf}
		\caption{Results of Singularity kernel with TSVD method and added white-noise.}
		\label{fig:Singularity_TSVD_Noise}
	\end{figure}
	
	\begin{figure}[h!]
		\centering
		\includegraphics[width=16cm]{../Outputs/Singularity_VCE_Noise.pdf}
		\caption{Results of Singularity kernel with Tikhonov (VCE) method and added white-noise.}
		\label{fig:Singularity_VCE_Noise}
	\end{figure}


	\subsection{Logarithmic With Noise}
	
	Results are shown in the figures \ref{fig:Logarithmic_Chol_Noise} to \ref{fig:Logarithmic_VCE_Noise}.Also, mean and norm ($L_2$) of the three methods are shown in table \ref{tab:Logarithmic_Error_Noise}.
	
	\begin{table}[h!]
		\centering
		\caption{Error for interpolation using Logarithmic kernel and added white-noise (unit of values are $\frac{m^2}{s^2}$).}
		\vspace{0.3cm}
		\renewcommand{\arraystretch}{1.4}
		\begin{tabular}{c|c|c|c}
			\textbf{Method} & Cholesky & TSVD & Tikhonov (VCE) \\
			\hline 
			\textbf{Mean Error} & 7.5111 & 7.3757 & 74.2444 \\
			\hline 
			\textbf{Norm of Errors} & 5155.9012 & 12404.6648 & 8167.3250 \\
		\end{tabular}
		\label{tab:Logarithmic_Error_Noise}
	\end{table}
	
	\begin{figure}[h!]
		\centering
		\includegraphics[width=16cm]{../Outputs/Logarithmic_Cholesky_Noise.pdf}
		\caption{Results of Logarithmic kernel with Cholesky decomposition and added white-noise.}
		\label{fig:Logarithmic_Chol_Noise}
	\end{figure}
	
	\clearpage
	
	\begin{figure}[h!]
		\centering
		\includegraphics[width=16cm]{../Outputs/Logarithmic_TSVD_Noise.pdf}
		\caption{Results of Logarithmic kernel with TSVD method and added white-noise.}
		\label{fig:Logarithmic_TSVD_Noise}
	\end{figure}
	
	\begin{figure}[h!]
		\centering
		\includegraphics[width=16cm]{../Outputs/Logarithmic_VCE_Noise.pdf}
		\caption{Results of Logarithmic kernel with Tikhonov (VCE) method and added white-noise.}
		\label{fig:Logarithmic_VCE_Noise}
	\end{figure}
	# Spherical-Spline-Interpolation


