// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// getSeqSimMatCpp
NumericMatrix getSeqSimMatCpp(std::string seq1, std::string seq2, double match, double misMatch);
RcppExport SEXP _DIAlignR_getSeqSimMatCpp(SEXP seq1SEXP, SEXP seq2SEXP, SEXP matchSEXP, SEXP misMatchSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type seq1(seq1SEXP);
    Rcpp::traits::input_parameter< std::string >::type seq2(seq2SEXP);
    Rcpp::traits::input_parameter< double >::type match(matchSEXP);
    Rcpp::traits::input_parameter< double >::type misMatch(misMatchSEXP);
    rcpp_result_gen = Rcpp::wrap(getSeqSimMatCpp(seq1, seq2, match, misMatch));
    return rcpp_result_gen;
END_RCPP
}
// getChromSimMatCpp
NumericMatrix getChromSimMatCpp(Rcpp::List l1, Rcpp::List l2, std::string normalization, std::string simType, double cosAngleThresh, double dotProdThresh, int kerLen);
RcppExport SEXP _DIAlignR_getChromSimMatCpp(SEXP l1SEXP, SEXP l2SEXP, SEXP normalizationSEXP, SEXP simTypeSEXP, SEXP cosAngleThreshSEXP, SEXP dotProdThreshSEXP, SEXP kerLenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type l1(l1SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type l2(l2SEXP);
    Rcpp::traits::input_parameter< std::string >::type normalization(normalizationSEXP);
    Rcpp::traits::input_parameter< std::string >::type simType(simTypeSEXP);
    Rcpp::traits::input_parameter< double >::type cosAngleThresh(cosAngleThreshSEXP);
    Rcpp::traits::input_parameter< double >::type dotProdThresh(dotProdThreshSEXP);
    Rcpp::traits::input_parameter< int >::type kerLen(kerLenSEXP);
    rcpp_result_gen = Rcpp::wrap(getChromSimMatCpp(l1, l2, normalization, simType, cosAngleThresh, dotProdThresh, kerLen));
    return rcpp_result_gen;
END_RCPP
}
// getGlobalAlignMaskCpp
NumericMatrix getGlobalAlignMaskCpp(const std::vector<double>& tA, const std::vector<double>& tB, const std::vector<double>& tBp, int noBeef, bool hardConstrain);
RcppExport SEXP _DIAlignR_getGlobalAlignMaskCpp(SEXP tASEXP, SEXP tBSEXP, SEXP tBpSEXP, SEXP noBeefSEXP, SEXP hardConstrainSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type tA(tASEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type tB(tBSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type tBp(tBpSEXP);
    Rcpp::traits::input_parameter< int >::type noBeef(noBeefSEXP);
    Rcpp::traits::input_parameter< bool >::type hardConstrain(hardConstrainSEXP);
    rcpp_result_gen = Rcpp::wrap(getGlobalAlignMaskCpp(tA, tB, tBp, noBeef, hardConstrain));
    return rcpp_result_gen;
END_RCPP
}
// constrainSimCpp
NumericMatrix constrainSimCpp(const NumericMatrix& sim, const NumericMatrix& MASK, double samples4gradient);
RcppExport SEXP _DIAlignR_constrainSimCpp(SEXP simSEXP, SEXP MASKSEXP, SEXP samples4gradientSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type sim(simSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type MASK(MASKSEXP);
    Rcpp::traits::input_parameter< double >::type samples4gradient(samples4gradientSEXP);
    rcpp_result_gen = Rcpp::wrap(constrainSimCpp(sim, MASK, samples4gradient));
    return rcpp_result_gen;
END_RCPP
}
// getBaseGapPenaltyCpp
double getBaseGapPenaltyCpp(const NumericMatrix& sim, std::string SimType, double gapQuantile);
RcppExport SEXP _DIAlignR_getBaseGapPenaltyCpp(SEXP simSEXP, SEXP SimTypeSEXP, SEXP gapQuantileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type sim(simSEXP);
    Rcpp::traits::input_parameter< std::string >::type SimType(SimTypeSEXP);
    Rcpp::traits::input_parameter< double >::type gapQuantile(gapQuantileSEXP);
    rcpp_result_gen = Rcpp::wrap(getBaseGapPenaltyCpp(sim, SimType, gapQuantile));
    return rcpp_result_gen;
END_RCPP
}
// areaIntegrator
NumericVector areaIntegrator(Rcpp::List l1, Rcpp::List l2, double left, double right, std::string integrationType, std::string baselineType, bool fitEMG, bool baseSubtraction, int kernelLen, int polyOrd);
RcppExport SEXP _DIAlignR_areaIntegrator(SEXP l1SEXP, SEXP l2SEXP, SEXP leftSEXP, SEXP rightSEXP, SEXP integrationTypeSEXP, SEXP baselineTypeSEXP, SEXP fitEMGSEXP, SEXP baseSubtractionSEXP, SEXP kernelLenSEXP, SEXP polyOrdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type l1(l1SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type l2(l2SEXP);
    Rcpp::traits::input_parameter< double >::type left(leftSEXP);
    Rcpp::traits::input_parameter< double >::type right(rightSEXP);
    Rcpp::traits::input_parameter< std::string >::type integrationType(integrationTypeSEXP);
    Rcpp::traits::input_parameter< std::string >::type baselineType(baselineTypeSEXP);
    Rcpp::traits::input_parameter< bool >::type fitEMG(fitEMGSEXP);
    Rcpp::traits::input_parameter< bool >::type baseSubtraction(baseSubtractionSEXP);
    Rcpp::traits::input_parameter< int >::type kernelLen(kernelLenSEXP);
    Rcpp::traits::input_parameter< int >::type polyOrd(polyOrdSEXP);
    rcpp_result_gen = Rcpp::wrap(areaIntegrator(l1, l2, left, right, integrationType, baselineType, fitEMG, baseSubtraction, kernelLen, polyOrd));
    return rcpp_result_gen;
END_RCPP
}
// sgolayCpp
NumericMatrix sgolayCpp(NumericMatrix chrom, int kernelLen, int polyOrd);
RcppExport SEXP _DIAlignR_sgolayCpp(SEXP chromSEXP, SEXP kernelLenSEXP, SEXP polyOrdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type chrom(chromSEXP);
    Rcpp::traits::input_parameter< int >::type kernelLen(kernelLenSEXP);
    Rcpp::traits::input_parameter< int >::type polyOrd(polyOrdSEXP);
    rcpp_result_gen = Rcpp::wrap(sgolayCpp(chrom, kernelLen, polyOrd));
    return rcpp_result_gen;
END_RCPP
}
// getAlignedTimesCpp
NumericMatrix getAlignedTimesCpp(Rcpp::List l1, Rcpp::List l2, int kernelLen, int polyOrd, std::string alignType, double adaptiveRT, std::string normalization, std::string simType, const std::vector<double>& Bp, double goFactor, double geFactor, double cosAngleThresh, bool OverlapAlignment, double dotProdThresh, double gapQuantile, int kerLen, bool hardConstrain, double samples4gradient);
RcppExport SEXP _DIAlignR_getAlignedTimesCpp(SEXP l1SEXP, SEXP l2SEXP, SEXP kernelLenSEXP, SEXP polyOrdSEXP, SEXP alignTypeSEXP, SEXP adaptiveRTSEXP, SEXP normalizationSEXP, SEXP simTypeSEXP, SEXP BpSEXP, SEXP goFactorSEXP, SEXP geFactorSEXP, SEXP cosAngleThreshSEXP, SEXP OverlapAlignmentSEXP, SEXP dotProdThreshSEXP, SEXP gapQuantileSEXP, SEXP kerLenSEXP, SEXP hardConstrainSEXP, SEXP samples4gradientSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type l1(l1SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type l2(l2SEXP);
    Rcpp::traits::input_parameter< int >::type kernelLen(kernelLenSEXP);
    Rcpp::traits::input_parameter< int >::type polyOrd(polyOrdSEXP);
    Rcpp::traits::input_parameter< std::string >::type alignType(alignTypeSEXP);
    Rcpp::traits::input_parameter< double >::type adaptiveRT(adaptiveRTSEXP);
    Rcpp::traits::input_parameter< std::string >::type normalization(normalizationSEXP);
    Rcpp::traits::input_parameter< std::string >::type simType(simTypeSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type Bp(BpSEXP);
    Rcpp::traits::input_parameter< double >::type goFactor(goFactorSEXP);
    Rcpp::traits::input_parameter< double >::type geFactor(geFactorSEXP);
    Rcpp::traits::input_parameter< double >::type cosAngleThresh(cosAngleThreshSEXP);
    Rcpp::traits::input_parameter< bool >::type OverlapAlignment(OverlapAlignmentSEXP);
    Rcpp::traits::input_parameter< double >::type dotProdThresh(dotProdThreshSEXP);
    Rcpp::traits::input_parameter< double >::type gapQuantile(gapQuantileSEXP);
    Rcpp::traits::input_parameter< int >::type kerLen(kerLenSEXP);
    Rcpp::traits::input_parameter< bool >::type hardConstrain(hardConstrainSEXP);
    Rcpp::traits::input_parameter< double >::type samples4gradient(samples4gradientSEXP);
    rcpp_result_gen = Rcpp::wrap(getAlignedTimesCpp(l1, l2, kernelLen, polyOrd, alignType, adaptiveRT, normalization, simType, Bp, goFactor, geFactor, cosAngleThresh, OverlapAlignment, dotProdThresh, gapQuantile, kerLen, hardConstrain, samples4gradient));
    return rcpp_result_gen;
END_RCPP
}
// alignChromatogramsCpp
S4 alignChromatogramsCpp(Rcpp::List l1, Rcpp::List l2, std::string alignType, const std::vector<double>& tA, const std::vector<double>& tB, std::string normalization, std::string simType, double B1p, double B2p, int noBeef, double goFactor, double geFactor, double cosAngleThresh, bool OverlapAlignment, double dotProdThresh, double gapQuantile, int kerLen, bool hardConstrain, double samples4gradient, std::string objType);
RcppExport SEXP _DIAlignR_alignChromatogramsCpp(SEXP l1SEXP, SEXP l2SEXP, SEXP alignTypeSEXP, SEXP tASEXP, SEXP tBSEXP, SEXP normalizationSEXP, SEXP simTypeSEXP, SEXP B1pSEXP, SEXP B2pSEXP, SEXP noBeefSEXP, SEXP goFactorSEXP, SEXP geFactorSEXP, SEXP cosAngleThreshSEXP, SEXP OverlapAlignmentSEXP, SEXP dotProdThreshSEXP, SEXP gapQuantileSEXP, SEXP kerLenSEXP, SEXP hardConstrainSEXP, SEXP samples4gradientSEXP, SEXP objTypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type l1(l1SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type l2(l2SEXP);
    Rcpp::traits::input_parameter< std::string >::type alignType(alignTypeSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type tA(tASEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type tB(tBSEXP);
    Rcpp::traits::input_parameter< std::string >::type normalization(normalizationSEXP);
    Rcpp::traits::input_parameter< std::string >::type simType(simTypeSEXP);
    Rcpp::traits::input_parameter< double >::type B1p(B1pSEXP);
    Rcpp::traits::input_parameter< double >::type B2p(B2pSEXP);
    Rcpp::traits::input_parameter< int >::type noBeef(noBeefSEXP);
    Rcpp::traits::input_parameter< double >::type goFactor(goFactorSEXP);
    Rcpp::traits::input_parameter< double >::type geFactor(geFactorSEXP);
    Rcpp::traits::input_parameter< double >::type cosAngleThresh(cosAngleThreshSEXP);
    Rcpp::traits::input_parameter< bool >::type OverlapAlignment(OverlapAlignmentSEXP);
    Rcpp::traits::input_parameter< double >::type dotProdThresh(dotProdThreshSEXP);
    Rcpp::traits::input_parameter< double >::type gapQuantile(gapQuantileSEXP);
    Rcpp::traits::input_parameter< int >::type kerLen(kerLenSEXP);
    Rcpp::traits::input_parameter< bool >::type hardConstrain(hardConstrainSEXP);
    Rcpp::traits::input_parameter< double >::type samples4gradient(samples4gradientSEXP);
    Rcpp::traits::input_parameter< std::string >::type objType(objTypeSEXP);
    rcpp_result_gen = Rcpp::wrap(alignChromatogramsCpp(l1, l2, alignType, tA, tB, normalization, simType, B1p, B2p, noBeef, goFactor, geFactor, cosAngleThresh, OverlapAlignment, dotProdThresh, gapQuantile, kerLen, hardConstrain, samples4gradient, objType));
    return rcpp_result_gen;
END_RCPP
}
// doAlignmentCpp
S4 doAlignmentCpp(NumericMatrix sim, double gap, bool OverlapAlignment);
RcppExport SEXP _DIAlignR_doAlignmentCpp(SEXP simSEXP, SEXP gapSEXP, SEXP OverlapAlignmentSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type sim(simSEXP);
    Rcpp::traits::input_parameter< double >::type gap(gapSEXP);
    Rcpp::traits::input_parameter< bool >::type OverlapAlignment(OverlapAlignmentSEXP);
    rcpp_result_gen = Rcpp::wrap(doAlignmentCpp(sim, gap, OverlapAlignment));
    return rcpp_result_gen;
END_RCPP
}
// doAffineAlignmentCpp
S4 doAffineAlignmentCpp(NumericMatrix sim, double go, double ge, bool OverlapAlignment);
RcppExport SEXP _DIAlignR_doAffineAlignmentCpp(SEXP simSEXP, SEXP goSEXP, SEXP geSEXP, SEXP OverlapAlignmentSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type sim(simSEXP);
    Rcpp::traits::input_parameter< double >::type go(goSEXP);
    Rcpp::traits::input_parameter< double >::type ge(geSEXP);
    Rcpp::traits::input_parameter< bool >::type OverlapAlignment(OverlapAlignmentSEXP);
    rcpp_result_gen = Rcpp::wrap(doAffineAlignmentCpp(sim, go, ge, OverlapAlignment));
    return rcpp_result_gen;
END_RCPP
}
// splineFillCpp
NumericVector splineFillCpp(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& xout);
RcppExport SEXP _DIAlignR_splineFillCpp(SEXP xSEXP, SEXP ySEXP, SEXP xoutSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type xout(xoutSEXP);
    rcpp_result_gen = Rcpp::wrap(splineFillCpp(x, y, xout));
    return rcpp_result_gen;
END_RCPP
}
// getChildXICpp
List getChildXICpp(Rcpp::List l1, Rcpp::List l2, int kernelLen, int polyOrd, std::string alignType, double adaptiveRT, std::string normalization, std::string simType, const std::vector<double>& Bp, double goFactor, double geFactor, double cosAngleThresh, bool OverlapAlignment, double dotProdThresh, double gapQuantile, int kerLen, bool hardConstrain, double samples4gradient, double wRef, std::string splineMethod, std::string mergeStrategy, bool keepFlanks);
RcppExport SEXP _DIAlignR_getChildXICpp(SEXP l1SEXP, SEXP l2SEXP, SEXP kernelLenSEXP, SEXP polyOrdSEXP, SEXP alignTypeSEXP, SEXP adaptiveRTSEXP, SEXP normalizationSEXP, SEXP simTypeSEXP, SEXP BpSEXP, SEXP goFactorSEXP, SEXP geFactorSEXP, SEXP cosAngleThreshSEXP, SEXP OverlapAlignmentSEXP, SEXP dotProdThreshSEXP, SEXP gapQuantileSEXP, SEXP kerLenSEXP, SEXP hardConstrainSEXP, SEXP samples4gradientSEXP, SEXP wRefSEXP, SEXP splineMethodSEXP, SEXP mergeStrategySEXP, SEXP keepFlanksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type l1(l1SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type l2(l2SEXP);
    Rcpp::traits::input_parameter< int >::type kernelLen(kernelLenSEXP);
    Rcpp::traits::input_parameter< int >::type polyOrd(polyOrdSEXP);
    Rcpp::traits::input_parameter< std::string >::type alignType(alignTypeSEXP);
    Rcpp::traits::input_parameter< double >::type adaptiveRT(adaptiveRTSEXP);
    Rcpp::traits::input_parameter< std::string >::type normalization(normalizationSEXP);
    Rcpp::traits::input_parameter< std::string >::type simType(simTypeSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type Bp(BpSEXP);
    Rcpp::traits::input_parameter< double >::type goFactor(goFactorSEXP);
    Rcpp::traits::input_parameter< double >::type geFactor(geFactorSEXP);
    Rcpp::traits::input_parameter< double >::type cosAngleThresh(cosAngleThreshSEXP);
    Rcpp::traits::input_parameter< bool >::type OverlapAlignment(OverlapAlignmentSEXP);
    Rcpp::traits::input_parameter< double >::type dotProdThresh(dotProdThreshSEXP);
    Rcpp::traits::input_parameter< double >::type gapQuantile(gapQuantileSEXP);
    Rcpp::traits::input_parameter< int >::type kerLen(kerLenSEXP);
    Rcpp::traits::input_parameter< bool >::type hardConstrain(hardConstrainSEXP);
    Rcpp::traits::input_parameter< double >::type samples4gradient(samples4gradientSEXP);
    Rcpp::traits::input_parameter< double >::type wRef(wRefSEXP);
    Rcpp::traits::input_parameter< std::string >::type splineMethod(splineMethodSEXP);
    Rcpp::traits::input_parameter< std::string >::type mergeStrategy(mergeStrategySEXP);
    Rcpp::traits::input_parameter< bool >::type keepFlanks(keepFlanksSEXP);
    rcpp_result_gen = Rcpp::wrap(getChildXICpp(l1, l2, kernelLen, polyOrd, alignType, adaptiveRT, normalization, simType, Bp, goFactor, geFactor, cosAngleThresh, OverlapAlignment, dotProdThresh, gapQuantile, kerLen, hardConstrain, samples4gradient, wRef, splineMethod, mergeStrategy, keepFlanks));
    return rcpp_result_gen;
END_RCPP
}
// otherChildXICpp
List otherChildXICpp(Rcpp::List l1, Rcpp::List l2, int kernelLen, int polyOrd, NumericMatrix mat, std::vector<double> childTime, double wRef, std::string splineMethod);
RcppExport SEXP _DIAlignR_otherChildXICpp(SEXP l1SEXP, SEXP l2SEXP, SEXP kernelLenSEXP, SEXP polyOrdSEXP, SEXP matSEXP, SEXP childTimeSEXP, SEXP wRefSEXP, SEXP splineMethodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type l1(l1SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type l2(l2SEXP);
    Rcpp::traits::input_parameter< int >::type kernelLen(kernelLenSEXP);
    Rcpp::traits::input_parameter< int >::type polyOrd(polyOrdSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type childTime(childTimeSEXP);
    Rcpp::traits::input_parameter< double >::type wRef(wRefSEXP);
    Rcpp::traits::input_parameter< std::string >::type splineMethod(splineMethodSEXP);
    rcpp_result_gen = Rcpp::wrap(otherChildXICpp(l1, l2, kernelLen, polyOrd, mat, childTime, wRef, splineMethod));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DIAlignR_getSeqSimMatCpp", (DL_FUNC) &_DIAlignR_getSeqSimMatCpp, 4},
    {"_DIAlignR_getChromSimMatCpp", (DL_FUNC) &_DIAlignR_getChromSimMatCpp, 7},
    {"_DIAlignR_getGlobalAlignMaskCpp", (DL_FUNC) &_DIAlignR_getGlobalAlignMaskCpp, 5},
    {"_DIAlignR_constrainSimCpp", (DL_FUNC) &_DIAlignR_constrainSimCpp, 3},
    {"_DIAlignR_getBaseGapPenaltyCpp", (DL_FUNC) &_DIAlignR_getBaseGapPenaltyCpp, 3},
    {"_DIAlignR_areaIntegrator", (DL_FUNC) &_DIAlignR_areaIntegrator, 10},
    {"_DIAlignR_sgolayCpp", (DL_FUNC) &_DIAlignR_sgolayCpp, 3},
    {"_DIAlignR_getAlignedTimesCpp", (DL_FUNC) &_DIAlignR_getAlignedTimesCpp, 18},
    {"_DIAlignR_alignChromatogramsCpp", (DL_FUNC) &_DIAlignR_alignChromatogramsCpp, 20},
    {"_DIAlignR_doAlignmentCpp", (DL_FUNC) &_DIAlignR_doAlignmentCpp, 3},
    {"_DIAlignR_doAffineAlignmentCpp", (DL_FUNC) &_DIAlignR_doAffineAlignmentCpp, 4},
    {"_DIAlignR_splineFillCpp", (DL_FUNC) &_DIAlignR_splineFillCpp, 3},
    {"_DIAlignR_getChildXICpp", (DL_FUNC) &_DIAlignR_getChildXICpp, 22},
    {"_DIAlignR_otherChildXICpp", (DL_FUNC) &_DIAlignR_otherChildXICpp, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_DIAlignR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
