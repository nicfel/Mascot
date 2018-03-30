#include "beast_mascot_ode_Euler2ndOrderNative.h"
#include "Euler2ndOrder.h"
#include "Euler2ndOrderCPU.h"
#include <iostream>

using namespace std;

//Euler2ndOrder * instance;


/*
 * Class:     beast_mascot_ode_Euler2ndOrderNative
 * Method:    setup
 * Signature: (I)V
 */
JNIEXPORT void JNICALL Java_beast_mascot_ode_Euler2ndOrderNative_setup
  (JNIEnv *env, jobject obj, jint maxSize){
//	instance = new Euler2ndOrderCPU();
//	instance->setup(maxSize);
	cerr << "Java_beast_mascot_ode_Euler2ndOrderNative_setup" << endl;
	cerr << "maxSize = " << maxSize << endl;
}

/*
 * Class:     beast_mascot_ode_Euler2ndOrderNative
 * Method:    init
 * Signature: ([[D[DIIDD)V
 */
JNIEXPORT void JNICALL Java_beast_mascot_ode_Euler2ndOrderNative_init
  (JNIEnv *env, jobject obj, jdoubleArray migration_ratesArray, jdoubleArray coalescent_ratesArray, jint lineages, jint states, jdouble epsilon, jdouble max_step){
//  	jdouble * migration_rates = (env)->GetDoubleArrayElements(migration_ratesArray, 0);
//  	jdouble * coalescent_rates = (env)->GetDoubleArrayElements(coalescent_ratesArray, 0);
//  	int rateCount = env->GetArrayLength(migration_ratesArray);
//
//	instance->init(migration_rates, rateCount, coalescent_rates, lineages, states, epsilon, max_step);
}

/*
 * Class:     beast_mascot_ode_Euler2ndOrderNative
 * Method:    initWithIndicators
 * Signature: ([[D[[I[DIIDD)V
 */
JNIEXPORT void JNICALL Java_beast_mascot_ode_Euler2ndOrderNative_initWithIndicators
  (JNIEnv *env, jobject obj, jdoubleArray, jobjectArray, jdoubleArray, jint, jint, jdouble, jdouble){
}

/*
 * Class:     beast_mascot_ode_Euler2ndOrderNative
 * Method:    calculateValues
 * Signature: (D[DI)V
 */
JNIEXPORT void JNICALL Java_beast_mascot_ode_Euler2ndOrderNative_calculateValues
  (JNIEnv *env, jobject obj, jdouble duration, jdoubleArray pArray, jint length){
//  	jdouble * p = (env)->GetDoubleArrayElements(pArray, 0);
//  	instance->calculateValues(duration, p, length);
//  	(env)->SetDoubleArrayRegion(pArray, 0, length, p);
}
