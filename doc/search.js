"use strict";
var items = [
{"atmosphere.mixture.MixtureOptimizerException" : "atmosphere/mixture/MixtureOptimizerException.html"},
{"atmosphere.mixture.MixtureOptimizerException.this" : "atmosphere/mixture/MixtureOptimizerException.this.html"},
{"atmosphere.mixture.FeaturesException" : "atmosphere/mixture/FeaturesException.html"},
{"atmosphere.mixture.FeaturesException.this" : "atmosphere/mixture/FeaturesException.this.html"},
{"atmosphere.mixture.MixtureOptimizer" : "atmosphere/mixture/MixtureOptimizer.html"},
{"atmosphere.mixture.MixtureOptimizer.this" : "atmosphere/mixture/MixtureOptimizer.this.html"},
{"atmosphere.mixture.MixtureOptimizer.eval" : "atmosphere/mixture/MixtureOptimizer.eval.html"},
{"atmosphere.mixture.MixtureOptimizer.update" : "atmosphere/mixture/MixtureOptimizer.update.html"},
{"atmosphere.mixture.MixtureOptimizer.optimize" : "atmosphere/mixture/MixtureOptimizer.optimize.html"},
{"atmosphere.mixture.MixtureOptimizer.optimize" : "atmosphere/mixture/MixtureOptimizer.optimize.html"},
{"atmosphere.mixture.MixtureOptimizer.optimize" : "atmosphere/mixture/MixtureOptimizer.optimize.html"},
{"atmosphere.mixture.MixtureOptimizer.put" : "atmosphere/mixture/MixtureOptimizer.put.html"},
{"atmosphere.mixture.MixtureOptimizer.put" : "atmosphere/mixture/MixtureOptimizer.put.html"},
{"atmosphere.mixture.MixtureOptimizer.length" : "atmosphere/mixture/MixtureOptimizer.length.html"},
{"atmosphere.mixture.MixtureOptimizer.maxLength" : "atmosphere/mixture/MixtureOptimizer.maxLength.html"},
{"atmosphere.mixture.MixtureOptimizer.reset" : "atmosphere/mixture/MixtureOptimizer.reset.html"},
{"atmosphere.mixture.MixtureOptimizer.popFront" : "atmosphere/mixture/MixtureOptimizer.popFront.html"},
{"atmosphere.mixture.MixtureOptimizer.popFrontN" : "atmosphere/mixture/MixtureOptimizer.popFrontN.html"},
{"atmosphere.mixture.MixtureOptimizer.features" : "atmosphere/mixture/MixtureOptimizer.features.html"},
{"atmosphere.mixture.MixtureOptimizer.mixture" : "atmosphere/mixture/MixtureOptimizer.mixture.html"},
{"atmosphere.mixture.MixtureOptimizer.weights" : "atmosphere/mixture/MixtureOptimizer.weights.html"},
{"atmosphere.mixture.MixtureOptimizer.weights" : "atmosphere/mixture/MixtureOptimizer.weights.html"},
{"atmosphere.mixture.GradientDescent" : "atmosphere/mixture/GradientDescent.html"},
{"atmosphere.mixture.GradientDescent.this" : "atmosphere/mixture/GradientDescent.this.html"},
{"atmosphere.mixture.CoordinateDescent" : "atmosphere/mixture/CoordinateDescent.html"},
{"atmosphere.mixture.CoordinateDescent.this" : "atmosphere/mixture/CoordinateDescent.this.html"},
{"atmosphere.mixture.CoordinateDescentPartial" : "atmosphere/mixture/CoordinateDescentPartial.html"},
{"atmosphere.mixture.CoordinateDescentPartial.this" : "atmosphere/mixture/CoordinateDescentPartial.this.html"},
{"atmosphere.mixture.LikelihoodMaximization" : "atmosphere/mixture/LikelihoodMaximization.html"},
{"atmosphere.mixture.LikelihoodMaximization.setWeightsInProportionToLikelihood" : "atmosphere/mixture/LikelihoodMaximization.setWeightsInProportionToLikelihood.html"},
{"atmosphere.mixture.LikelihoodMaximization.put" : "atmosphere/mixture/LikelihoodMaximization.put.html"},
{"atmosphere.mixture.LikelihoodMaximization.putAndSetWeightsInProportionToLikelihood" : "atmosphere/mixture/LikelihoodMaximization.putAndSetWeightsInProportionToLikelihood.html"},
{"atmosphere.mixture.LikelihoodMaximization.optimize" : "atmosphere/mixture/LikelihoodMaximization.optimize.html"},
{"atmosphere.mixture.LikelihoodMaximization.optimize" : "atmosphere/mixture/LikelihoodMaximization.optimize.html"},
{"atmosphere.mixture.LikelihoodMaximization.isFeaturesCorrect" : "atmosphere/mixture/LikelihoodMaximization.isFeaturesCorrect.html"},
{"atmosphere.mixture.LikelihoodMaximization.log2Likelihood" : "atmosphere/mixture/LikelihoodMaximization.log2Likelihood.html"},
{"atmosphere.mixture.CoordinateLikelihoodMaximization" : "atmosphere/mixture/CoordinateLikelihoodMaximization.html"},
{"atmosphere.mixture.CoordinateLikelihoodMaximization.this" : "atmosphere/mixture/CoordinateLikelihoodMaximization.this.html"},
{"atmosphere.mixture.GradientLikelihoodMaximization" : "atmosphere/mixture/GradientLikelihoodMaximization.html"},
{"atmosphere.mixture.GradientLikelihoodMaximization.this" : "atmosphere/mixture/GradientLikelihoodMaximization.this.html"},
{"atmosphere.parametrized.nvmm.NormalVarianceMeanMixture" : "atmosphere/parametrized/nvmm/NormalVarianceMeanMixture.html"},
{"atmosphere.parametrized.nvmm.NormalVarianceMeanMixture.this" : "atmosphere/parametrized/nvmm/NormalVarianceMeanMixture.this.html"},
{"atmosphere.parametrized.nvmm.NormalVarianceMeanMixture.optimize" : "atmosphere/parametrized/nvmm/NormalVarianceMeanMixture.optimize.html"},
{"atmosphere.parametrized.nvmm.NormalVarianceMeanMixture.optimize" : "atmosphere/parametrized/nvmm/NormalVarianceMeanMixture.optimize.html"},
{"atmosphere.parametrized.nvmm.NormalVarianceMeanMixture.optimize" : "atmosphere/parametrized/nvmm/NormalVarianceMeanMixture.optimize.html"},
{"atmosphere.parametrized.nvmm.NormalVarianceMeanMixture.sample" : "atmosphere/parametrized/nvmm/NormalVarianceMeanMixture.sample.html"},
{"atmosphere.parametrized.nvmm.NormalVarianceMeanMixture.sample" : "atmosphere/parametrized/nvmm/NormalVarianceMeanMixture.sample.html"},
{"atmosphere.parametrized.nvmm.NormalVarianceMeanMixture.mean" : "atmosphere/parametrized/nvmm/NormalVarianceMeanMixture.mean.html"},
{"atmosphere.parametrized.nvmm.NormalVarianceMeanMixture.alpha" : "atmosphere/parametrized/nvmm/NormalVarianceMeanMixture.alpha.html"},
{"atmosphere.parametrized.nvmm.NormalVarianceMeanMixture.grid" : "atmosphere/parametrized/nvmm/NormalVarianceMeanMixture.grid.html"},
{"atmosphere.parametrized.nvmm.NormalVarianceMeanMixtureEM" : "atmosphere/parametrized/nvmm/NormalVarianceMeanMixtureEM.html"},
{"atmosphere.parametrized.nvmm.NormalVarianceMeanMixtureEM.this" : "atmosphere/parametrized/nvmm/NormalVarianceMeanMixtureEM.this.html"},
{"atmosphere.parametrized.nvmm.NormalVarianceMeanMixtureEMAndGradient" : "atmosphere/parametrized/nvmm/NormalVarianceMeanMixtureEMAndGradient.html"},
{"atmosphere.parametrized.nvmm.NormalVarianceMeanMixtureEMAndGradient.this" : "atmosphere/parametrized/nvmm/NormalVarianceMeanMixtureEMAndGradient.this.html"},
{"atmosphere.parametrized.nvmm.NormalVarianceMeanMixtureEMAndCoordinate" : "atmosphere/parametrized/nvmm/NormalVarianceMeanMixtureEMAndCoordinate.html"},
{"atmosphere.parametrized.nvmm.NormalVarianceMeanMixtureEMAndCoordinate.this" : "atmosphere/parametrized/nvmm/NormalVarianceMeanMixtureEMAndCoordinate.this.html"},
{"distribution.cdf.CDF" : "distribution/cdf/CDF.html"},
{"distribution.cdf.CDF.opCall" : "distribution/cdf/CDF.opCall.html"},
{"distribution.cdf.GammaCDF" : "distribution/cdf/GammaCDF.html"},
{"distribution.cdf.GammaCDF.this" : "distribution/cdf/GammaCDF.this.html"},
{"distribution.cdf.InverseGammaCDF" : "distribution/cdf/InverseGammaCDF.html"},
{"distribution.cdf.InverseGammaCDF.this" : "distribution/cdf/InverseGammaCDF.this.html"},
{"distribution.pdf.PDF" : "distribution/pdf/PDF.html"},
{"distribution.pdf.PDF.opCall" : "distribution/pdf/PDF.opCall.html"},
{"distribution.pdf.GammaPDF" : "distribution/pdf/GammaPDF.html"},
{"distribution.pdf.GammaPDF.this" : "distribution/pdf/GammaPDF.this.html"},
{"distribution.pdf.InverseGammaPDF" : "distribution/pdf/InverseGammaPDF.html"},
{"distribution.pdf.InverseGammaPDF.this" : "distribution/pdf/InverseGammaPDF.this.html"},
{"distribution.pdf.GeneralizedGammaPDF" : "distribution/pdf/GeneralizedGammaPDF.html"},
{"distribution.pdf.GeneralizedGammaPDF.this" : "distribution/pdf/GeneralizedGammaPDF.this.html"},
{"distribution.pdf.InverseGaussianPDF" : "distribution/pdf/InverseGaussianPDF.html"},
{"distribution.pdf.InverseGaussianPDF.this" : "distribution/pdf/InverseGaussianPDF.this.html"},
{"distribution.pdf.ProperGeneralizedInverseGaussianPDF" : "distribution/pdf/ProperGeneralizedInverseGaussianPDF.html"},
{"distribution.pdf.ProperGeneralizedInverseGaussianPDF.this" : "distribution/pdf/ProperGeneralizedInverseGaussianPDF.this.html"},
{"distribution.pdf.GeneralizedInverseGaussianPDF" : "distribution/pdf/GeneralizedInverseGaussianPDF.html"},
{"distribution.pdf.GeneralizedInverseGaussianPDF.this" : "distribution/pdf/GeneralizedInverseGaussianPDF.this.html"},
{"distribution.pdf.VarianceGammaPDF" : "distribution/pdf/VarianceGammaPDF.html"},
{"distribution.pdf.VarianceGammaPDF.this" : "distribution/pdf/VarianceGammaPDF.this.html"},
{"distribution.pdf.HyperbolicAsymmetricTPDF" : "distribution/pdf/HyperbolicAsymmetricTPDF.html"},
{"distribution.pdf.HyperbolicAsymmetricTPDF.this" : "distribution/pdf/HyperbolicAsymmetricTPDF.this.html"},
{"distribution.pdf.NormalInverseGaussianPDF" : "distribution/pdf/NormalInverseGaussianPDF.html"},
{"distribution.pdf.NormalInverseGaussianPDF.this" : "distribution/pdf/NormalInverseGaussianPDF.this.html"},
{"distribution.pdf.ProperGeneralizedHyperbolicPDF" : "distribution/pdf/ProperGeneralizedHyperbolicPDF.html"},
{"distribution.pdf.ProperGeneralizedHyperbolicPDF.this" : "distribution/pdf/ProperGeneralizedHyperbolicPDF.this.html"},
{"distribution.pdf.GeneralizedHyperbolicPDF" : "distribution/pdf/GeneralizedHyperbolicPDF.html"},
{"distribution.pdf.GeneralizedHyperbolicPDF.this" : "distribution/pdf/GeneralizedHyperbolicPDF.this.html"},
{"distribution.quantile.Quantile" : "distribution/quantile/Quantile.html"},
{"distribution.quantile.Quantile.opCall" : "distribution/quantile/Quantile.opCall.html"},
{"distribution.quantile.GammaQuantile" : "distribution/quantile/GammaQuantile.html"},
{"distribution.quantile.GammaQuantile.this" : "distribution/quantile/GammaQuantile.this.html"},
{"distribution.quantile.InverseGammaQuantile" : "distribution/quantile/InverseGammaQuantile.html"},
{"distribution.quantile.InverseGammaQuantile.this" : "distribution/quantile/InverseGammaQuantile.this.html"},
{"distribution.random.DistributionRNG" : "distribution/random/DistributionRNG.html"},
{"distribution.random.DistributionRNG.empty" : "distribution/random/DistributionRNG.empty.html"},
{"distribution.random.DistributionRNG.popFront" : "distribution/random/DistributionRNG.popFront.html"},
{"distribution.random.DistributionRNG.front" : "distribution/random/DistributionRNG.front.html"},
{"distribution.random.NormalVarianceMeanMixtureRNG" : "distribution/random/NormalVarianceMeanMixtureRNG.html"},
{"distribution.random.NormalVarianceMeanMixtureRNG.this" : "distribution/random/NormalVarianceMeanMixtureRNG.this.html"},
{"distribution.random.GammaRNG" : "distribution/random/GammaRNG.html"},
{"distribution.random.GammaRNG.this" : "distribution/random/GammaRNG.this.html"},
{"distribution.random.InverseGammaRNG" : "distribution/random/InverseGammaRNG.html"},
{"distribution.random.InverseGammaRNG.this" : "distribution/random/InverseGammaRNG.this.html"},
{"distribution.random.GeneralizedGammaRNG" : "distribution/random/GeneralizedGammaRNG.html"},
{"distribution.random.GeneralizedGammaRNG.this" : "distribution/random/GeneralizedGammaRNG.this.html"},
{"distribution.random.InverseGaussianRNG" : "distribution/random/InverseGaussianRNG.html"},
{"distribution.random.InverseGaussianRNG.this" : "distribution/random/InverseGaussianRNG.this.html"},
{"distribution.random.ProperGeneralizedInverseGaussianRNG" : "distribution/random/ProperGeneralizedInverseGaussianRNG.html"},
{"distribution.random.ProperGeneralizedInverseGaussianRNG.this" : "distribution/random/ProperGeneralizedInverseGaussianRNG.this.html"},
{"distribution.random.GeneralizedInverseGaussianRNG" : "distribution/random/GeneralizedInverseGaussianRNG.html"},
{"distribution.random.GeneralizedInverseGaussianRNG.this" : "distribution/random/GeneralizedInverseGaussianRNG.this.html"},
{"distribution.random.VarianceGammaRNG" : "distribution/random/VarianceGammaRNG.html"},
{"distribution.random.VarianceGammaRNG.this" : "distribution/random/VarianceGammaRNG.this.html"},
{"distribution.random.HyperbolicAsymmetricTRNG" : "distribution/random/HyperbolicAsymmetricTRNG.html"},
{"distribution.random.HyperbolicAsymmetricTRNG.this" : "distribution/random/HyperbolicAsymmetricTRNG.this.html"},
{"distribution.random.GeneralizedVarianceGammaRNG" : "distribution/random/GeneralizedVarianceGammaRNG.html"},
{"distribution.random.GeneralizedVarianceGammaRNG.this" : "distribution/random/GeneralizedVarianceGammaRNG.this.html"},
{"distribution.random.NormalInverseGaussianRNG" : "distribution/random/NormalInverseGaussianRNG.html"},
{"distribution.random.NormalInverseGaussianRNG.this" : "distribution/random/NormalInverseGaussianRNG.this.html"},
{"distribution.random.ProperGeneralizedHyperbolicRNG" : "distribution/random/ProperGeneralizedHyperbolicRNG.html"},
{"distribution.random.ProperGeneralizedHyperbolicRNG.this" : "distribution/random/ProperGeneralizedHyperbolicRNG.this.html"},
{"distribution.random.GeneralizedHyperbolicRNG" : "distribution/random/GeneralizedHyperbolicRNG.html"},
{"distribution.random.GeneralizedHyperbolicRNG.this" : "distribution/random/GeneralizedHyperbolicRNG.this.html"},
{"distribution.random.uniformM11E0" : "distribution/random/uniformM11E0.html"},
{"distribution.random.uniformM11E0" : "distribution/random/uniformM11E0.html"},
{"distribution.random.rNormal" : "distribution/random/rNormal.html"},
{"distribution.random.rNormal" : "distribution/random/rNormal.html"},
{"distribution.random.rExponential" : "distribution/random/rExponential.html"},
{"distribution.random.rExponential" : "distribution/random/rExponential.html"},
{"distribution.random.rGamma" : "distribution/random/rGamma.html"},
{"distribution.random.rGamma" : "distribution/random/rGamma.html"},
{"distribution.random.rGeneralizedGamma" : "distribution/random/rGeneralizedGamma.html"},
{"distribution.random.rGeneralizedGamma" : "distribution/random/rGeneralizedGamma.html"},
{"distribution.random.rInverseGamma" : "distribution/random/rInverseGamma.html"},
{"distribution.random.rInverseGamma" : "distribution/random/rInverseGamma.html"},
{"distribution.random.rInverseGaussian" : "distribution/random/rInverseGaussian.html"},
{"distribution.random.rInverseGaussian" : "distribution/random/rInverseGaussian.html"},
{"distribution.random.rChiSquare" : "distribution/random/rChiSquare.html"},
{"distribution.random.rChiSquare" : "distribution/random/rChiSquare.html"},
{"distribution.random.rStudentT" : "distribution/random/rStudentT.html"},
{"distribution.random.rStudentT" : "distribution/random/rStudentT.html"},
{"distribution.random.rWeibull" : "distribution/random/rWeibull.html"},
{"distribution.random.rWeibull" : "distribution/random/rWeibull.html"},
{"simple_matrix.Vector" : "simple_matrix/Vector.html"},
{"simple_matrix.Vector.ptr" : "simple_matrix/Vector.ptr.html"},
{"simple_matrix.Vector.length" : "simple_matrix/Vector.length.html"},
{"simple_matrix.Vector.shift" : "simple_matrix/Vector.shift.html"},
{"simple_matrix.Vector.popFront" : "simple_matrix/Vector.popFront.html"},
{"simple_matrix.Vector.popFrontN" : "simple_matrix/Vector.popFrontN.html"},
{"simple_matrix.Vector.back" : "simple_matrix/Vector.back.html"},
{"simple_matrix.Vector.popBack" : "simple_matrix/Vector.popBack.html"},
{"simple_matrix.Vector.popBackN" : "simple_matrix/Vector.popBackN.html"},
{"simple_matrix.Vector.opIndex" : "simple_matrix/Vector.opIndex.html"},
{"simple_matrix.Vector.opIndex" : "simple_matrix/Vector.opIndex.html"},
{"simple_matrix.Matrix" : "simple_matrix/Matrix.html"},
{"simple_matrix.Matrix.ptr" : "simple_matrix/Matrix.ptr.html"},
{"simple_matrix.Matrix.height" : "simple_matrix/Matrix.height.html"},
{"simple_matrix.Matrix.width" : "simple_matrix/Matrix.width.html"},
{"simple_matrix.Matrix.shift" : "simple_matrix/Matrix.shift.html"},
{"simple_matrix.Matrix.ptrEnd" : "simple_matrix/Matrix.ptrEnd.html"},
{"simple_matrix.Matrix.transposed" : "simple_matrix/Matrix.transposed.html"},
{"simple_matrix.Matrix.transversal" : "simple_matrix/Matrix.transversal.html"},
{"simple_matrix.Matrix.lengthTransveral" : "simple_matrix/Matrix.lengthTransveral.html"},
{"simple_matrix.Matrix.emptyTransveral" : "simple_matrix/Matrix.emptyTransveral.html"},
{"simple_matrix.Matrix.frontTransversal" : "simple_matrix/Matrix.frontTransversal.html"},
{"simple_matrix.Matrix.popFrontTransversal" : "simple_matrix/Matrix.popFrontTransversal.html"},
{"simple_matrix.Matrix.popFrontNTransversal" : "simple_matrix/Matrix.popFrontNTransversal.html"},
{"simple_matrix.Matrix.backTransversal" : "simple_matrix/Matrix.backTransversal.html"},
{"simple_matrix.Matrix.popBackTransversal" : "simple_matrix/Matrix.popBackTransversal.html"},
{"simple_matrix.Matrix.popBackNTransversal" : "simple_matrix/Matrix.popBackNTransversal.html"},
{"simple_matrix.Matrix.transpose" : "simple_matrix/Matrix.transpose.html"},
{"simple_matrix.Matrix.this" : "simple_matrix/Matrix.this.html"},
{"simple_matrix.Matrix.this" : "simple_matrix/Matrix.this.html"},
{"simple_matrix.Matrix.this" : "simple_matrix/Matrix.this.html"},
{"simple_matrix.Matrix.this" : "simple_matrix/Matrix.this.html"},
{"simple_matrix.Matrix.front" : "simple_matrix/Matrix.front.html"},
{"simple_matrix.Matrix.popFront" : "simple_matrix/Matrix.popFront.html"},
{"simple_matrix.Matrix.popFrontN" : "simple_matrix/Matrix.popFrontN.html"},
{"simple_matrix.Matrix.back" : "simple_matrix/Matrix.back.html"},
{"simple_matrix.Matrix.popBack" : "simple_matrix/Matrix.popBack.html"},
{"simple_matrix.Matrix.popBackN" : "simple_matrix/Matrix.popBackN.html"},
{"simple_matrix.Matrix.opIndex" : "simple_matrix/Matrix.opIndex.html"},
{"simple_matrix.Matrix.opIndex" : "simple_matrix/Matrix.opIndex.html"},
{"simple_matrix.Matrix.opIndex" : "simple_matrix/Matrix.opIndex.html"},
{"simple_matrix.TransposedMatrix" : "simple_matrix/TransposedMatrix.html"},
{"simple_matrix.TransposedMatrix.matrix" : "simple_matrix/TransposedMatrix.matrix.html"},
{"simple_matrix.TransposedMatrix.width" : "simple_matrix/TransposedMatrix.width.html"},
{"simple_matrix.TransposedMatrix.height" : "simple_matrix/TransposedMatrix.height.html"},
{"simple_matrix.TransposedMatrix.lengthTransveral" : "simple_matrix/TransposedMatrix.lengthTransveral.html"},
{"simple_matrix.TransposedMatrix.emptyTransveral" : "simple_matrix/TransposedMatrix.emptyTransveral.html"},
{"simple_matrix.TransposedMatrix.transposed" : "simple_matrix/TransposedMatrix.transposed.html"},
{"simple_matrix.TransposedMatrix.frontTransversal" : "simple_matrix/TransposedMatrix.frontTransversal.html"},
{"simple_matrix.TransposedMatrix.popFrontTransversal" : "simple_matrix/TransposedMatrix.popFrontTransversal.html"},
{"simple_matrix.TransposedMatrix.popFrontNTransversal" : "simple_matrix/TransposedMatrix.popFrontNTransversal.html"},
{"simple_matrix.TransposedMatrix.backTransversal" : "simple_matrix/TransposedMatrix.backTransversal.html"},
{"simple_matrix.TransposedMatrix.popBackTransversal" : "simple_matrix/TransposedMatrix.popBackTransversal.html"},
{"simple_matrix.TransposedMatrix.popBackNTransversal" : "simple_matrix/TransposedMatrix.popBackNTransversal.html"},
{"simple_matrix.TransposedMatrix.front" : "simple_matrix/TransposedMatrix.front.html"},
{"simple_matrix.TransposedMatrix.popFront" : "simple_matrix/TransposedMatrix.popFront.html"},
{"simple_matrix.TransposedMatrix.popFrontN" : "simple_matrix/TransposedMatrix.popFrontN.html"},
{"simple_matrix.TransposedMatrix.back" : "simple_matrix/TransposedMatrix.back.html"},
{"simple_matrix.TransposedMatrix.popBack" : "simple_matrix/TransposedMatrix.popBack.html"},
{"simple_matrix.TransposedMatrix.popBackN" : "simple_matrix/TransposedMatrix.popBackN.html"},
{"simple_matrix.TransposedMatrix.opIndex" : "simple_matrix/TransposedMatrix.opIndex.html"},
{"simple_matrix.TransposedMatrix.opIndex" : "simple_matrix/TransposedMatrix.opIndex.html"},
{"simple_matrix.TransposedMatrix.opIndex" : "simple_matrix/TransposedMatrix.opIndex.html"},
{"simple_matrix.TransposedMatrix.opIndex" : "simple_matrix/TransposedMatrix.opIndex.html"},
{"simple_matrix.SlidingWindow" : "simple_matrix/SlidingWindow.html"},
{"simple_matrix.SlidingWindow.data" : "simple_matrix/SlidingWindow.data.html"},
{"simple_matrix.SlidingWindow.transposedMatrix" : "simple_matrix/SlidingWindow.transposedMatrix.html"},
{"simple_matrix.SlidingWindow.this" : "simple_matrix/SlidingWindow.this.html"},
{"simple_matrix.SlidingWindow.reset" : "simple_matrix/SlidingWindow.reset.html"},
{"simple_matrix.SlidingWindow.put" : "simple_matrix/SlidingWindow.put.html"},
{"simple_matrix.SlidingWindow.reserveBackN" : "simple_matrix/SlidingWindow.reserveBackN.html"},
{"simple_matrix.SlidingWindow.moveToFront" : "simple_matrix/SlidingWindow.moveToFront.html"},
{"simple_matrix._D1" : "simple_matrix/_D1.html"},
{"simple_matrix._D1.empty" : "simple_matrix/_D1.empty.html"},
{"simple_matrix._D1.save" : "simple_matrix/_D1.save.html"},
{"simple_matrix._D1.opIndex" : "simple_matrix/_D1.opIndex.html"},
{"simple_matrix._D1.opDollar" : "simple_matrix/_D1.opDollar.html"},
{"simple_matrix._D1.opSlice" : "simple_matrix/_D1.opSlice.html"},
{"simple_matrix._D2" : "simple_matrix/_D2.html"},
{"simple_matrix._D2.length" : "simple_matrix/_D2.length.html"},
{"simple_matrix._D2.opDollar" : "simple_matrix/_D2.opDollar.html"},
{"simple_matrix._D2.opDollar" : "simple_matrix/_D2.opDollar.html"},
];
function search(str) {
	var re = new RegExp(str.toLowerCase());
	var ret = {};
	for (var i = 0; i < items.length; i++) {
		var k = Object.keys(items[i])[0];
		if (re.test(k.toLowerCase()))
			ret[k] = items[i][k];
	}
	return ret;
}

function searchSubmit(value, event) {
	console.log("searchSubmit");
	var resultTable = document.getElementById("results");
	while (resultTable.firstChild)
		resultTable.removeChild(resultTable.firstChild);
	if (value === "" || event.keyCode == 27) {
		resultTable.style.display = "none";
		return;
	}
	resultTable.style.display = "block";
	var results = search(value);
	var keys = Object.keys(results);
	if (keys.length === 0) {
		var row = resultTable.insertRow();
		var td = document.createElement("td");
		var node = document.createTextNode("No results");
		td.appendChild(node);
		row.appendChild(td);
		return;
	}
	for (var i = 0; i < keys.length; i++) {
		var k = keys[i];
		var v = results[keys[i]];
		var link = document.createElement("a");
		link.href = v;
		link.textContent = k;
		link.attributes.id = "link" + i;
		var row = resultTable.insertRow();
		row.appendChild(link);
	}
}

function hideSearchResults(event) {
	if (event.keyCode != 27)
		return;
	var resultTable = document.getElementById("results");
	while (resultTable.firstChild)
		resultTable.removeChild(resultTable.firstChild);
	resultTable.style.display = "none";
}

