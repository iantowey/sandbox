package com.lucidlogic.horsetracker.handler;

import android.view.View;

import com.lucidlogic.horsetracker.model.Race;
import com.lucidlogic.horsetracker.model.Racecard;

import timber.log.Timber;

/**
 * Created by ian on 12/07/16.
 */
public class RacecardClickHandler {

    public void onRacecardClick(View v, Racecard r){
        Timber.i("%s clicked" ,r.getDate());
    }
}
