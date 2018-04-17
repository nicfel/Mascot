#include "beast_mascot_distribution_MascotNative2.h"
#include "Mascot.h"
#include <stdio.h>

Mascot  * instance;

/*
 * Class:     beast_mascot_distribution_MascotNative2
 * Method:    setup
 * Signature: ([IIDDIII)V
 */
JNIEXPORT void JNICALL Java_beast_mascot_distribution_MascotNative2_setup
  (JNIEnv * env, jobject o, jintArray nodeTypeArray, jint states,
		  jdouble epsilon, jdouble max_step, jint sampleCount, jint nodeCount, jint intervalCount, jboolean useCache) {
	jint * nodeType = (env)->GetIntArrayElements(nodeTypeArray, 0);
	instance = new Mascot(nodeType, states, epsilon, max_step, sampleCount, nodeCount, nodeCount, useCache);
}

/*
 * Class:     beast_mascot_distribution_MascotNative2
 * Method:    calculateLogP
 * Signature: (ZI[I[I[D[I)D
 */
JNIEXPORT jdouble JNICALL Java_beast_mascot_distribution_MascotNative2_calculateLogP
  (JNIEnv * env, jobject o, jboolean dynamicsIsDirty, jint firstDirtyInterval,
		  jintArray lineagesAddedArray, jintArray lineagesRemovedArray, jdoubleArray intervalsArray,
		  jintArray parentsArray) {
	jint * lineagesAdded = (env)->GetIntArrayElements(lineagesAddedArray, 0);
	jint * lineagesRemoved = (env)->GetIntArrayElements(lineagesRemovedArray, 0);
	jdouble * intervals = (env)->GetDoubleArrayElements(intervalsArray, 0);
	jint * parents = (env)->GetIntArrayElements(parentsArray, 0);

	jdouble logP = instance->calculateLogP(dynamicsIsDirty, firstDirtyInterval,
			lineagesAdded, lineagesRemoved, intervals, parents);
	//fprintf(stderr, "logP = %f\n", logP);
	return logP;
}

/*
 * Class:     beast_mascot_distribution_MascotNative2
 * Method:    setUpDynamics
 * Signature: ([D[D[[I[D)V
 */
JNIEXPORT void JNICALL Java_beast_mascot_distribution_MascotNative2_setUpDynamics
  (JNIEnv * env, jobject o, jdoubleArray coalescentRatesArray, jdoubleArray migrationRatesArray, jobjectArray indicatorsMatrix, jdoubleArray nextRateShiftArray) {
	jdouble * coalescentRates = (env)->GetDoubleArrayElements(coalescentRatesArray, 0);
	jdouble * migrationRates = (env)->GetDoubleArrayElements(migrationRatesArray, 0);
	jdouble * nextRateShift = (env)->GetDoubleArrayElements(nextRateShiftArray, 0);
	int rateShiftCount = env->GetArrayLength(nextRateShiftArray);

//	  	fprintf(stderr, "ratesInterval %d\n", rateShiftCount);
//	  	for (int k = 0;k < rateShiftCount; k++) {
//	  		fprintf(stderr, "%d %f\n", k, nextRateShift[k]);
//	  	for (int i = 0; i < 5; i++) {
//	  		fprintf(stderr, "[");
//	  		for (int j = 0; j < 5; j++) {
//	  			fprintf(stderr, "%f ", migrationRates[k * 25 + i*5 + j]);
//	  		}
//	  		fprintf(stderr, "] %f \n", coalescentRates[k * 5 + i]);
//	  	}
//	  	}

	instance->setUpDynamics(rateShiftCount, coalescentRates, migrationRates, nextRateShift);
}

/*
 * Class:     beast_mascot_distribution_MascotNative2
 * Method:    getStateProb
 * Signature: (I[D)V
 */
JNIEXPORT void JNICALL Java_beast_mascot_distribution_MascotNative2_getStateProb
  (JNIEnv * env, jobject o, jint nr, jdoubleArray pArray) {
	jdouble * p = (env)->GetDoubleArrayElements(pArray, 0);
	instance->getStateProb(nr, p);
  	(env)->SetDoubleArrayRegion(pArray, 0, instance->states, p);
}

/*
 * Class:     beast_mascot_distribution_MascotNative2
 * Method:    storeState
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_beast_mascot_distribution_MascotNative2_storeState
  (JNIEnv * env, jobject o) {
	instance->store();
}

/*
 * Class:     beast_mascot_distribution_MascotNative2
 * Method:    restoreState
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_beast_mascot_distribution_MascotNative2_restoreState
  (JNIEnv * env, jobject o) {
	instance->restore();
}
