#include "beast_mascot_ode_Euler2ndOrderNative.h"
#include "Euler2ndOrder.h"
#include "Euler2ndOrderCPU.h"
#include <stdio.h>

Euler2ndOrderCPU * instance;


/*
 * Class:     beast_mascot_ode_Euler2ndOrderNative
 * Method:    setup
 * Signature: (I)V
 */
JNIEXPORT void JNICALL Java_beast_mascot_ode_Euler2ndOrderNative_setup
  (JNIEnv *env, jobject obj, jint maxSize, jint states, jdouble epsilon, jdouble max_step) {
	switch (states) {
	case 2: instance = new Euler2ndOrderCPU2();break;
	case 3: instance = new Euler2ndOrderCPU3();break;
	case 4: instance = new Euler2ndOrderCPU4();break;
	case 5: instance = new Euler2ndOrderCPU5();break;
	case 6: instance = new Euler2ndOrderCPU6();break;
	case 7: instance = new Euler2ndOrderCPU7();break;
	case 8: instance = new Euler2ndOrderCPU8();break;
	case 9: instance = new Euler2ndOrderCPU9();break;
	case 10: instance = new Euler2ndOrderCPU10();break;
	default: instance = new Euler2ndOrderCPU();
	}
	instance->setup(maxSize, states, epsilon, max_step);
	printf("Java_beast_mascot_ode_Euler2ndOrderNative_setup\n");
	printf("maxSize = %d\n", maxSize);
}

/*
 * Class:     beast_mascot_ode_Euler2ndOrderNative
 * Method:    init
 * Signature: ([[D[DIIDD)V
 */
JNIEXPORT void JNICALL Java_beast_mascot_ode_Euler2ndOrderNative_init
  (JNIEnv *env, jobject obj, jdoubleArray migration_ratesArray, jdoubleArray coalescent_ratesArray, jint lineages) {
  	jdouble * migration_rates = (env)->GetDoubleArrayElements(migration_ratesArray, 0);
  	jdouble * coalescent_rates = (env)->GetDoubleArrayElements(coalescent_ratesArray, 0);
  	int rateCount = env->GetArrayLength(migration_ratesArray);

	instance->init(migration_rates, rateCount, coalescent_rates, lineages);
}

/*
 * Class:     beast_mascot_ode_Euler2ndOrderNative
 * Method:    initWithIndicators
 * Signature: ([[D[[I[DIIDD)V
 */
JNIEXPORT void JNICALL Java_beast_mascot_ode_Euler2ndOrderNative_initWithIndicators
  (JNIEnv *env, jobject obj, jdoubleArray migration_ratesArray, jintArray indicatorsArray, jdoubleArray coalescent_ratesArray, jint lineages) {
}

/*
 * Class:     beast_mascot_ode_Euler2ndOrderNative
 * Method:    calculateValues
 * Signature: (D[DI)V
 */
JNIEXPORT void JNICALL Java_beast_mascot_ode_Euler2ndOrderNative_calculateValues
  (JNIEnv *env, jobject obj, jdouble duration, jdoubleArray pArray, jint length) {
  	jdouble * p = (env)->GetDoubleArrayElements(pArray, 0);
  	instance->calculateValues(duration, p, length);
  	(env)->SetDoubleArrayRegion(pArray, 0, length, p);
}


/*
 * Class:     beast_mascot_ode_Euler2ndOrderNative
 * Method:    initAndcalculateValues
 * Signature: ([D[DIIDDD[DI)V
 */
JNIEXPORT void JNICALL Java_beast_mascot_ode_Euler2ndOrderNative_initAndcalculateValues
(JNIEnv *env, jobject obj, jdoubleArray migration_ratesArray, jdoubleArray coalescent_ratesArray, jint lineages,
		jdouble duration, jdoubleArray pArray, jint length) {
  	jdouble * migration_rates = (env)->GetDoubleArrayElements(migration_ratesArray, 0);
  	jdouble * coalescent_rates = (env)->GetDoubleArrayElements(coalescent_ratesArray, 0);
  	int rateCount = env->GetArrayLength(migration_ratesArray);

	instance->init(migration_rates, rateCount, coalescent_rates, lineages);

	jdouble * p = (env)->GetDoubleArrayElements(pArray, 0);
  	instance->calculateValues(duration, p, length);
  	(env)->SetDoubleArrayRegion(pArray, 0, length, p);
}
