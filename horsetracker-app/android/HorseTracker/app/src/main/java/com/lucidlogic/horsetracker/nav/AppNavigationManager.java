package com.lucidlogic.horsetracker.nav;

import android.support.v4.app.Fragment;
import android.support.v4.app.FragmentManager;
import android.support.v7.app.AppCompatActivity;

import com.lucidlogic.horsetracker.R;

import timber.log.Timber;

/**
 * Created by ian on 01/08/16.
 */
public class AppNavigationManager {

    private static AppCompatActivity _activity;

    public static void initialize(AppCompatActivity activity){
        _activity = activity;
    }

    public static void addFragment(Fragment fragment, String label) {
        _activity.getSupportFragmentManager()
                .beginTransaction()
                .replace(R.id.main_display, fragment, label)
                .addToBackStack(label)
                .commit();
        Timber.i("fragment with label %s added" , label);
    }

    public static Fragment findFragment(String tag) {
        Fragment frag = _activity.getSupportFragmentManager().findFragmentByTag(tag);
        return frag;
    }

    public static FragmentManager getFragmentManager(){
        return _activity.getSupportFragmentManager();
    }
}
