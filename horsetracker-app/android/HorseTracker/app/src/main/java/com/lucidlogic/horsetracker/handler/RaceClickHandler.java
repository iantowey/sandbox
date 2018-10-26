package com.lucidlogic.horsetracker.handler;

import android.support.v4.app.Fragment;
import android.view.View;

import com.lucidlogic.horsetracker.R;
import com.lucidlogic.horsetracker.config.Constants;
import com.lucidlogic.horsetracker.model.Entity;
import com.lucidlogic.horsetracker.model.Race;
import com.lucidlogic.horsetracker.nav.AppNavigationManager;
import com.lucidlogic.horsetracker.view.impl.RaceFragment;

import timber.log.Timber;

/**
 * Created by ian on 12/07/16.
 */
public class RaceClickHandler {

    public void onRaceClick(View v, Race r){
        Timber.i("%s %s clicked" ,r.getTime(), r.getName());
        String fragmentTag = Constants.RACE + r.getId();
        Fragment fragment = AppNavigationManager.findFragment(fragmentTag);
        if(fragment == null){
            fragment = RaceFragment.newInstance(r.getId());
        }
        AppNavigationManager.addFragment(fragment, fragmentTag);
    }
}