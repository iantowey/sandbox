package com.lucidlogic.horsetracker.handler;

import android.support.v4.app.Fragment;
import android.view.View;

import com.lucidlogic.horsetracker.config.Constants;
import com.lucidlogic.horsetracker.model.Entity;
import com.lucidlogic.horsetracker.nav.AppNavigationManager;
import com.lucidlogic.horsetracker.view.impl.EntityFragment;

import timber.log.Timber;

/**
 * Created by ian on 12/07/16.
 */
public class RunnerEntityClickHandler {

    public void onHorseClick(View v, Entity e){
        Timber.i("%s %s clicked" ,e.getName(), e.getId());
        String fragmentTag = Constants.ENTITY + e.getId();
        Fragment fragment = AppNavigationManager.findFragment(fragmentTag);
        if(fragment == null){
            fragment = EntityFragment.newInstance(e.getId());
        }
        AppNavigationManager.addFragment(fragment, fragmentTag);
    }

    public void onJockeyClick(View v, Entity entity){
        Timber.i("JOCKEY :: %s %s clicked" ,entity.getName(), entity.getId());
    }

}
