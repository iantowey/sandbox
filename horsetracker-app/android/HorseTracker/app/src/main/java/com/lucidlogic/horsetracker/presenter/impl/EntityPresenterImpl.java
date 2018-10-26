package com.lucidlogic.horsetracker.presenter.impl;

import com.lucidlogic.horsetracker.presenter.EntityPresenter;
import com.lucidlogic.horsetracker.service.ParseService;
import com.lucidlogic.horsetracker.utils.BeanTransformers;
import com.lucidlogic.horsetracker.view.EntityView;

import rx.android.schedulers.AndroidSchedulers;
import rx.schedulers.Schedulers;
import timber.log.Timber;

/**
 * Created by itowey on 11/08/16.
 */
public class EntityPresenterImpl implements EntityPresenter {

    private  EntityView entityView;

    @Override
    public void setView(EntityView entityView) {
        this.entityView =  entityView;
    }

    @Override
    public void updateView(String entityObjectId) {
        ParseService.getEntity(entityObjectId)
                .observeOn(AndroidSchedulers.mainThread())
                .subscribeOn(Schedulers.io())
                .map(entityDTO -> BeanTransformers.entityFromEntityDTO(entityDTO, true))
                .subscribe(entity -> {
                            Timber.i("Item retrieved");
                            entityView.updateView(entity);
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
