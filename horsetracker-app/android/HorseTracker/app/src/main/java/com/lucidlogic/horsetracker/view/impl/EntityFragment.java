package com.lucidlogic.horsetracker.view.impl;


import android.databinding.DataBindingUtil;
import android.databinding.ViewDataBinding;
import android.os.Bundle;
import android.support.v4.app.Fragment;
import android.view.LayoutInflater;
import android.view.View;
import android.view.ViewGroup;
import android.widget.Button;

import com.annimon.stream.Stream;
import com.lucidlogic.horsetracker.BR;
import com.lucidlogic.horsetracker.R;
import com.lucidlogic.horsetracker.config.Constants;
import com.lucidlogic.horsetracker.handler.EntityClickHandler;
import com.lucidlogic.horsetracker.model.Entity;
import com.lucidlogic.horsetracker.model.dto.UserFollowsEntityDTO;
import com.lucidlogic.horsetracker.presenter.EntityPresenter;
import com.lucidlogic.horsetracker.presenter.impl.EntityPresenterImpl;
import com.lucidlogic.horsetracker.service.ParseService;
import com.lucidlogic.horsetracker.view.EntityView;

import java.util.List;

/**
 * A simple {@link Fragment} subclass.
 */
public class EntityFragment extends Fragment implements EntityView {

    public static final String ENTITY_OBJECT_ID = "ENTITY_OBJECT_ID";
    private String entityObjectId;
    private EntityPresenter entityPresenter;
    private Entity entity;
    private EntityClickHandler entityClickHandler;

    public EntityFragment() {
        entityPresenter = new EntityPresenterImpl();
        entityPresenter.setView(this);
        this.entity = new Entity();
        this.entityClickHandler = new EntityClickHandler();
    }

    public static EntityFragment newInstance(String entityObjectId) {
        EntityFragment fragment = new EntityFragment();
        Bundle args = new Bundle();
        args.putString(ENTITY_OBJECT_ID, entityObjectId);
        fragment.setArguments(args);
        return fragment;
    }

    @Override
    public void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        if (getArguments() != null) {
            entityObjectId = getArguments().getString(ENTITY_OBJECT_ID);
        }
    }


    @Override
    public View onCreateView(LayoutInflater inflater, ViewGroup container, Bundle savedInstanceState) {

        ViewDataBinding viewDataBinding = DataBindingUtil.inflate(inflater, R.layout.fragment_entity, container, false);

        viewDataBinding.setVariable(BR.entity,  entity);
        viewDataBinding.setVariable(BR.handler,  new EntityClickHandler());
        entityPresenter.updateView(this.entityObjectId);
        View v = viewDataBinding.getRoot();
        Button b = (Button)v.findViewById(R.id.entityFollowBtn);
        b.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                entityClickHandler.onEntityClick(entity);
            }
        });
        return v;
    }

    @Override
    public void updateView(Entity e) {
        this.entity = e;
        List<UserFollowsEntityDTO>  userFollowsEntityDTOs = ParseService.getFromLocalParseStoreByLabel(Constants.USER_FOLLOWS_ENTITY_LIST, UserFollowsEntityDTO.class);
        boolean userTracking = Stream.of(userFollowsEntityDTOs).anyMatch(value -> value.getEntity().getObjectId().equals(entity.getId()));
        entity.setUserFollowing(userTracking);
        this.entity.notifyChange();
    }
}
