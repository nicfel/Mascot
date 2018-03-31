/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class beast_mascot_ode_Euler2ndOrderNative */

#ifndef _Included_beast_mascot_ode_Euler2ndOrderNative
#define _Included_beast_mascot_ode_Euler2ndOrderNative
#ifdef __cplusplus
extern "C" {
#endif
/*
 * Class:     beast_mascot_ode_Euler2ndOrderNative
 * Method:    setup
 * Signature: (I)V
 */
JNIEXPORT void JNICALL Java_beast_mascot_ode_Euler2ndOrderNative_setup
  (JNIEnv *, jobject, jint);

/*
 * Class:     beast_mascot_ode_Euler2ndOrderNative
 * Method:    init
 * Signature: ([D[DIIDD)V
 */
JNIEXPORT void JNICALL Java_beast_mascot_ode_Euler2ndOrderNative_init
  (JNIEnv *, jobject, jdoubleArray, jdoubleArray, jint, jint, jdouble, jdouble);

/*
 * Class:     beast_mascot_ode_Euler2ndOrderNative
 * Method:    initWithIndicators
 * Signature: ([D[I[DIIDD)V
 */
JNIEXPORT void JNICALL Java_beast_mascot_ode_Euler2ndOrderNative_initWithIndicators
  (JNIEnv *, jobject, jdoubleArray, jintArray, jdoubleArray, jint, jint, jdouble, jdouble);

/*
 * Class:     beast_mascot_ode_Euler2ndOrderNative
 * Method:    calculateValues
 * Signature: (D[DI)V
 */
JNIEXPORT void JNICALL Java_beast_mascot_ode_Euler2ndOrderNative_calculateValues
  (JNIEnv *, jobject, jdouble, jdoubleArray, jint);

/*
 * Class:     beast_mascot_ode_Euler2ndOrderNative
 * Method:    initAndcalculateValues
 * Signature: ([D[DIIDDD[DI)V
 */
JNIEXPORT void JNICALL Java_beast_mascot_ode_Euler2ndOrderNative_initAndcalculateValues
  (JNIEnv *, jobject, jdoubleArray, jdoubleArray, jint, jint, jdouble, jdouble, jdouble, jdoubleArray, jint);

#ifdef __cplusplus
}
#endif
#endif