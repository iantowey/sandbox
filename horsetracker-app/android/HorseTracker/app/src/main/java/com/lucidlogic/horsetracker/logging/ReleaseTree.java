package com.lucidlogic.horsetracker.logging;


import android.util.Log;
import com.crashlytics.android.Crashlytics;
import timber.log.Timber;

public class ReleaseTree extends Timber.Tree {

    private static final String CRASHLYTICS_KEY_PRIORITY = "priority";
    private static final String CRASHLYTICS_KEY_TAG = "tag";
    private static final String CRASHLYTICS_KEY_MESSAGE = "message";

    @Override
    protected boolean isLoggable(int priority) {
        return !(priority == Log.VERBOSE || priority == Log.DEBUG || priority == Log.INFO);
    }

    @Override
    protected void log(int priority, String tag, String message, Throwable t) {
        if (isLoggable(priority)) { //Not sure if this is needed...
            Timber.i("logging exception....");

            Crashlytics.setInt(CRASHLYTICS_KEY_PRIORITY, priority);
            Crashlytics.setString(CRASHLYTICS_KEY_TAG, tag);
            Crashlytics.setString(CRASHLYTICS_KEY_MESSAGE, message);

            if (t == null) {
                Crashlytics.logException(new Exception(message));
            } else {
                Crashlytics.logException(t);
            }
        }
    }
}