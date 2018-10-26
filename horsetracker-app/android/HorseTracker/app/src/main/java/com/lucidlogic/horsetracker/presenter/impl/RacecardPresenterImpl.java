package com.lucidlogic.horsetracker.presenter.impl;

import com.lucidlogic.horsetracker.model.Racecard;
import com.lucidlogic.horsetracker.model.dto.RacecardDTO;
import com.lucidlogic.horsetracker.presenter.RacecardPresenter;
import com.lucidlogic.horsetracker.service.ParseService;
import com.lucidlogic.horsetracker.utils.BeanTransformers;
import com.lucidlogic.horsetracker.view.RacecardView;

import java.util.Date;

import rx.Observer;
import rx.android.schedulers.AndroidSchedulers;
import rx.functions.Func1;
import rx.schedulers.Schedulers;
import timber.log.Timber;

/**
 * Created by ian on 10/07/16.
 */
public class RacecardPresenterImpl implements RacecardPresenter {
    RacecardView racecardView;

    @Override
    public void setView(RacecardView racecardView) {

        this.racecardView = racecardView;
    }

    @Override
    public void updateView() {
        ParseService.getRacecards(new Date())
                .observeOn(AndroidSchedulers.mainThread())
                .subscribeOn(Schedulers.io())
                .map(racecardDTO -> BeanTransformers.racecardFromRacecardDTO(racecardDTO, false))
                .subscribe(racecard -> {
                            Timber.i("Item retrieved");
                            racecardView.updateView(racecard);
                        },
                        e -> {
                            Timber.e(e,e.getMessage());
                            e.printStackTrace();
                        },
                        () -> {
                            Timber.i("Item retrival complete");
                        });
    }
}
