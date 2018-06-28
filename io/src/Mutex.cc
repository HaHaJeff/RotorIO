#include "Mutex.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static void PthreadCall(const char* label, int result) {
	if (result != 0) {
		fprintf(stderr, "pthread %s: %s\n", label, strerror(result));
		// If the abort() function casues process termination, all open streams are closed and flushed;
		abort();
	}
}

Mutex::Mutex() { PthreadCall("init mutex", pthread_mutex_init(&mu_, NULL)); }

Mutex::~Mutex() { PthreadCall("destroy mutex", pthread_mutex_destroy(&mu_)); }

void Mutex::Lock() { PthreadCall("lock", pthread_mutex_lock(&mu_)); }

void Mutex::Unlock() { PthreadCall("unlock", pthread_mutex_unlock(&mu_)); }

CondVar::CondVar(Mutex* mu)
	: mu_(mu) {
		PthreadCall("init cv", pthread_cond_init(&cv_, NULL));
}

CondVar::~CondVar() { PthreadCall("destroy cv", pthread_cond_destroy(&cv_)); }

void CondVar::Wait() {
	PthreadCall("wait", pthread_cond_wait(&cv_, &mu_->mu_));
}

void CondVar::Signal() {
	PthreadCall("signal", pthread_cond_signal(&cv_));
}

void CondVar::SignalAll() {
	PthreadCall("broadcasr", pthread_cond_broadcast(&cv_));
}

void InitOnce(OnceType* once, void (*initializer)()) {
	PthreadCall("once", pthread_once(once, initializer));
}
